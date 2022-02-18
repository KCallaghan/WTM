//#define EIGEN_DONT_PARALLELIZE

#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.
#include <eigen3/Eigen/Core>
#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"
#include "mat_mult.cpp"
#include "add_recharge.hpp"


using namespace Eigen;

typedef Eigen::SparseMatrix<double,RowMajor> SpMat; // declares a row-major sparse matrix type of double
typedef Eigen::Triplet<double> T;  // used to populate the matrices

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {


double depthIntegratedTransmissivity(
  const double wtd_T,
  const double fdepth,
  const double ksat
){
  constexpr double shallow = 1.5;
  //Global soil datasets include information for shallow soils.
  //if the water table is deeper than this, the permeability
  //of the soil sees an exponential decay with depth.
  if(fdepth<=0) {
    std::cout<<"if"<<std::endl;
    // If the fdepth is zero, there is no water transmission below the surface
    // soil layer.
    // If it is less than zero, it is incorrect -- but no water transmission
    // also seems an okay thing to do in this case.
    return 0;
  } else if(wtd_T < -shallow){ // Equation S6 from the Fan paper
    return std::max(0.0,fdepth * ksat  * std::exp((wtd_T + shallow)/fdepth));
  } else if(wtd_T > 0){
    // If wtd_T is greater than 0, max out rate of groundwater movement
    // as though wtd_T were 0. The surface water will get to move in
    // FillSpillMerge.
    return std::max(0.0,ksat * (0 + shallow + fdepth));
  } else { //Equation S4 from the Fan paper
    return std::max(0.0,ksat * (wtd_T + shallow + fdepth));  //max because you can't have a negative transmissivity.
  }
}


//update the entire transmissivity array
void updateTransmissivity(
  Parameters &params,ArrayPack &arp,int continue_picard
){

  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) != 0.f){
        arp.my_last_wtd(x,y) = (arp.wtd_T(x,y) + arp.wtd_T_iteration(x,y))/2.;
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.my_last_wtd(x,y), arp.fdepth(x,y), arp.ksat(x,y));
      }
    }
  }




int first_half(const Parameters &params,ArrayPack &arp, int picard_number){

  std::vector<T> coefficients;
  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x*params.ncells_y);

  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);


  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);

    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
      //populate the known vector b. This is the current wtd. We also populate 'guess', which is based on the previous picard iteration's wtd result and we will use it as the guess to get our answer, x.
      b(y+(x*params.ncells_y)) = 0.;
      guess(y+(x*params.ncells_y)) = 0.;
    }
    else{
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
      guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);
    }
  }

  double entry;
  int main_loc = 0;

  //HALFWAY SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=1;x<params.ncells_x; x++)
  for(int y=1;y<params.ncells_y; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain, provided that they are not ocean cells.
//    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
//      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
//      coefficients.push_back(T(main_loc,main_loc, entry));
//    }
 //   else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(arp.scalar_array_y(x,y)/(2*arp.effective_storativity(x,y))) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(params.x_partial/(2*arp.effective_storativity(x,y))) +1;
      coefficients.push_back(T(main_loc,main_loc, entry));

      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      if(y != params.ncells_y-1 && y!= 0){  // && arp.land_mask(x,y+1) != 0.f){
        //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
        entry = -arp.scalar_array_y(x,y)/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x,y+1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+1,entry ));
      }
      //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
      if(y != 0 && y != params.ncells_y-1){  // && arp.land_ mask(x,y-1) != 0.f){
        //y may not == 0 since then there is no cell to the west.
        entry = -arp.scalar_array_y(x,y)/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x,y-1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-1,entry ));
      }
      //Now let's do the North diagonal. Offset by -(ncells_y).
      if(x != 0 && x != params.ncells_x-1){  // && arp.land_mask(x-1,y) != 0.f){
        //x may not equal 0 since then there is no cell to the north.
        entry = -params.x_partial/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x-1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-params.ncells_y, entry));
      }
      //finally, do the South diagonal, offset by +(ncells_y).
      if(x != params.ncells_x-1 && x!=0){  // && arp.land_mask(x+1,y) != 0.f){
        //we may not be in the final row where there is no cell
        entry = -params.x_partial/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x+1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+params.ncells_y, entry));
      }
  //  }
  }

  //use the triplet vector to populate the matrix, A.
  A.setFromTriplets(coefficients.begin(),coefficients.end());


  // Biconjugate gradient solver with guess
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(0.00001);
  //NOTE: we cannot use the Eigen:IncompleteLUT preconditioner, because its implementation is serial. Using it means that BiCGSTAB will not run in parallel. It is faster without.

  std::cout<<"compute"<<std::endl;
  solver.compute(A);

  assert(solver.info()==Eigen::Success);

  std::cout<<"solve"<<std::endl;
  vec_x = solver.solveWithGuess(b, guess);

  assert(solver.info()==Eigen::Success);

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  int cell_count = 0;

  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){

  //copy result into the wtd_T array:
    if(arp.land_mask(x,y) != 0.f){
      arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
    }
    else
      arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - 0.;


  }

  if(picard_number == 3)
    return 0;
  else{
    picard_number += 1;
    return picard_number;
  }
 // }

}






void second_half(Parameters &params,ArrayPack &arp){

  std::vector<T> coefficients_A;
  std::vector<T> coefficients_B;
  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);
  SpMat B(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);

  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x*params.ncells_y);


  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){

    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
      //populate the known vector b. This is the current wtd. We also populate 'guess', which is based on the previous picard iteration's wtd result and we will use it as the guess to get our answer, x.
      b(y+(x*params.ncells_y)) = 0.;
      guess(y+(x*params.ncells_y)) = 0.;
    }
    else{
      b(y+(x*params.ncells_y)) = arp.original_wtd(x,y) + arp.topo(x,y);
      guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);
    }
  }

  double entry;
  int main_loc = 0;


  ////SECOND SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=1;x<params.ncells_x; x++)
  for(int y=1;y<params.ncells_y; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain, provided that they are not ocean cells.
//    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
//      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
//      coefficients_A.push_back(T(main_loc,main_loc, entry));
//      coefficients_B.push_back(T(main_loc,main_loc, entry));
//    }
//    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      entry =  1 - (arp.scalar_array_x(x,y)*(arp.transmissivity(x-1,y)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x+1,y)/4.) ) - (arp.scalar_array_y(x,y)*(arp.transmissivity(x,y-1)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x,y+1)/4.));
      coefficients_B.push_back(T(main_loc,main_loc, entry));
      entry =  1 + (arp.scalar_array_x(x,y)*(arp.transmissivity(x-1,y)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x+1,y)/4.)) + (arp.scalar_array_y(x,y)*(arp.transmissivity(x,y-1)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x,y+1)/4.));
      coefficients_A.push_back(T(main_loc,main_loc, entry));

     //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      if(y != params.ncells_y-1 && y!= 0){  // && arp.land_mask(x,y+1) != 0.f){
        //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
        entry = (arp.scalar_array_y(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x,y+1)/4.);
        coefficients_B.push_back(T(main_loc,main_loc+1,entry ));
        entry = (-arp.scalar_array_y(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x,y+1)/4.);
        coefficients_A.push_back(T(main_loc,main_loc+1,entry ));
      }
      //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
      if(y != 0 && y != params.ncells_y-1){  // && arp.land_ mask(x,y-1) != 0.f){
        //y may not == 0 since then there is no cell to the west.
        entry = (arp.scalar_array_y(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x,y-1)/4.);
        coefficients_B.push_back(T(main_loc,main_loc-1,entry ));
        entry = (-arp.scalar_array_y(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x,y-1)/4.);
        coefficients_A.push_back(T(main_loc,main_loc-1,entry ));
      }
      //Now let's do the North diagonal. Offset by -(ncells_y).
      if(x != 0 && x != params.ncells_x-1){  // && arp.land_mask(x-1,y) != 0.f){
        //x may not equal 0 since then there is no cell to the north.
        entry = (arp.scalar_array_x(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x-1,y)/4.);
        coefficients_B.push_back(T(main_loc,main_loc-params.ncells_y,entry ));
        entry = (-arp.scalar_array_x(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x-1,y)/4.);
        coefficients_A.push_back(T(main_loc,main_loc-params.ncells_y,entry ));
      }
      //finally, do the South diagonal, offset by +(ncells_y).
      if(x != params.ncells_x-1 && x!=0){  // && arp.land_mask(x+1,y) != 0.f){
        //we may not be in the final row where there is no cell
        entry = (arp.scalar_array_x(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x+1,y)/4.);
        coefficients_B.push_back(T(main_loc,main_loc+params.ncells_y,entry ));
        entry = (-arp.scalar_array_x(x,y))*(arp.transmissivity(x,y)/4. + arp.transmissivity(x+1,y)/4.);
        coefficients_A.push_back(T(main_loc,main_loc+params.ncells_y,entry ));
        }
  //  }
  }

  std::cerr<<"finished second set of matrices"<<std::endl;

  //use the triplet vector to populate the matrices, A and B.
  A.setFromTriplets(coefficients_A.begin(),coefficients_A.end());
  B.setFromTriplets(coefficients_B.begin(),coefficients_B.end());


  b = B*b;
  guess = B*guess;


 // Apply the second half of the recharge to the water-table depth grid (wtd_T is currently in use). For this, start with the original wtd and add the second half of rech to it.
  //TODO: is this right, or should we start with the halfway-solved wtd_T and add half the rech to it?
 //   #pragma omp parallel for collapse(2)
    for(int y=0;y<params.ncells_y;y++)
    for(int x=0;x<params.ncells_x; x++){
      double rech_change = arp.rech(x,y)/31536000. * params.deltat;
      if(arp.original_wtd(x,y) >= 0 && arp.land_mask(x,y) == 1){          //skip ocean cells
        params.total_added_recharge += rech_change*arp.cell_area[y];
          b(y+(x*params.ncells_y)) += rech_change;
          guess(y+(x*params.ncells_y)) += rech_change;
          if(arp.original_wtd(x,y) + rech_change <0){
            double temp = -(arp.original_wtd(x,y) + rech_change);
            b(y+(x*params.ncells_y)) += temp;
            guess(y+(x*params.ncells_y)) += temp;
            params.total_added_recharge += temp*arp.cell_area[y];  //in this scenario, there has been a negative amount of recharge, so temp here will be positive and remove the extra amount that was spuriously subtracted.
          }
        }
        else if (rech_change>0 && arp.land_mask(x,y) == 1){
          params.total_added_recharge += rech_change*arp.cell_area[y];
          double GW_space = -arp.original_wtd(x,y) * arp.porosity(x,y);
          if(GW_space > rech_change){
            b(y+(x*params.ncells_y)) += rech_change/arp.porosity(x,y);
            guess(y+(x*params.ncells_y)) += rech_change/arp.porosity(x,y);
          }
          else{
            double temp = ( rech_change - GW_space) - arp.original_wtd(x,y);
            b(y+(x*params.ncells_y)) += temp;
            guess(y+(x*params.ncells_y)) += temp;
          }
        }
      }





  // Biconjugate gradient solver with guess
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(0.00001);
  //NOTE: we cannot use the Eigen:IncompleteLUT preconditioner, because its implementation is serial. Using it means that BiCGSTAB will not run in parallel. It is faster without.

  std::cout<<"compute"<<std::endl;
  solver.compute(A);

  assert(solver.info()==Eigen::Success);

  std::cout<<"solve"<<std::endl;
  vec_x = solver.solveWithGuess(b, guess);

  assert(solver.info()==Eigen::Success);

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;


  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
  //copy result into the wtd_T array:
    if(arp.land_mask(x,y) != 0.f){
      arp.wtd(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);

    }
    else
      arp.wtd(x,y) = vec_x(y+(x*params.ncells_y)) - 0.;


  }
  std::cerr<<"finished assigning the new wtd"<<std::endl;
}





//Deal with the fact that porosity is 1 above ground but [value] below ground. It is included directly in the matrix calculation,
//so we can't just scale water change before or after. Instead, we are scaling the effective storativity according to how much
//change in water table we are projecting to see during that time step.
void updateEffectiveStorativity(const Parameters &params,ArrayPack &arp){
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    float projected_full_step_wtd = arp.original_wtd(x,y) + (arp.wtd_T(x,y) - arp.original_wtd(x,y))*2.;
    if(arp.original_wtd(x,y) <= 0.f && projected_full_step_wtd <= 0.f) //both are below ground, so we can use the original porosity
      arp.effective_storativity(x,y) = arp.porosity(x,y);
    else if(arp.original_wtd(x,y) >= 0.f && projected_full_step_wtd >= 0.f) //both are above ground, so the porosity is 1
      arp.effective_storativity(x,y) = 1.;
    else{
      float change_in_water_column_thickness = std::abs(projected_full_step_wtd - arp.original_wtd(x,y));
      if(arp.original_wtd(x,y) < 0.f && projected_full_step_wtd > 0.f){ //started below ground and ended above ground
        // First, scale the change in wtd as if the whole column has the porosity of the belowground area
        float scaled_change_in_water_column_thickness = change_in_water_column_thickness * arp.effective_storativity(x,y)/arp.porosity(x,y);
        // Then get just that portion that is >0 (aboveground), and scale it back down to the equivalent thickness with 100% porosity (surface water)
        float aboveground_water_column_thickness = (scaled_change_in_water_column_thickness + arp.original_wtd(x,y)) * arp.effective_storativity(x,y);//Above-ground effective porosity is ==1, so no need to actually divide by 1 here (save on computation).
        //belowground water column thickness is equal to -original_wtd, so no need to assign it to a new variable.
        //float belowground_water_column_thickness = arp.original_wtd(x,y);//this is what the other formula works out to, I think. No more space than this below the ground. //scaled_change_in_water_column_thickness - (scaled_change_in_water_column_thickness + arp.original_wtd(x,y));
        arp.effective_storativity(x,y) = (aboveground_water_column_thickness + arp.porosity(x,y)*-arp.original_wtd(x,y)) / (aboveground_water_column_thickness - arp.original_wtd(x,y));
      }
      else if(arp.original_wtd(x,y) > 0.f && projected_full_step_wtd < 0.f){ //started above ground and ended below ground
        float scaled_change_in_water_column_thickness = change_in_water_column_thickness * arp.effective_storativity(x,y); //divided by 1 for above-ground porosity
        if(scaled_change_in_water_column_thickness < projected_full_step_wtd){  //when rescaling the water according to the porosity values, there is no longer enough to reach below ground; so the scaled new porosity will be = 1
          arp.effective_storativity(x,y) = 1;
        }
        else{
          // Get the belowground water thickness and expand it to a deeper depth in order to account for changing porosity
          float belowground_water_column_thickness = (scaled_change_in_water_column_thickness - arp.original_wtd(x,y)) * arp.effective_storativity(x,y)/arp.porosity(x,y);
          arp.effective_storativity(x,y) = (arp.original_wtd(x,y) + arp.porosity(x,y)*belowground_water_column_thickness) / (arp.original_wtd(x,y) + belowground_water_column_thickness);
        }
      }
    }
  }
}










//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void UpdateCPU(Parameters &params, ArrayPack &arp){

  Eigen::initParallel();
  omp_set_num_threads(8);
  Eigen::setNbThreads(8);

  // Picard iteration through solver
  params.x_partial = params.deltat/(params.cellsize_n_s_metres*params.cellsize_n_s_metres);
  float ocean_T = 0.000001 * (1.5 + 25); //some constant for all T values in ocean - TODO look up representative values

std::cout<<"create some needed arrays "<<std::endl;

// Apply the first half of the recharge to the water-table depth grid (wtd)
    // Its clone (wtd_T) is used and updated in the Picard iteration
    //also set the starting porosity
//set the scalar arrays for x and y directions
  //#pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      arp.transmissivity(x,y) = ocean_T;
      arp.temp_T(x,y) = ocean_T;
      params.total_loss_to_ocean += arp.wtd(x,y)*arp.cell_area[y];
      arp.wtd(x,y) = 0.;
      arp.wtd_T(x,y) = 0.;
      arp.original_wtd(x,y) = 0.;
      arp.wtd_T_iteration(x,y) = 0.;
      arp.effective_storativity(x,y) = 1.;
      }
    else{
      arp.original_wtd(x,y) = arp.wtd(x,y);
      arp.wtd(x,y) = add_recharge(params.deltat/2., arp.rech(x,y), arp.wtd(x,y), arp.land_mask(x,y), arp.porosity(x,y)); //use regular porosity for adding recharge since this checks for underground space within add_recharge.
      arp.wtd_T(x,y) = arp.wtd(x,y);
      arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);

      if(arp.original_wtd(x,y)>0)
        arp.effective_storativity(x,y) = 1;
      else
        arp.effective_storativity(x,y) = arp.porosity(x,y);

    }
      arp.scalar_array_y(x,y) = params.deltat/(arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]);
  }



   int continue_picard = 1;
   while(continue_picard != 0 ){
    std::cout << "updateTransmissivity: " << std::endl;
    updateTransmissivity(params,arp,continue_picard);
    std::cout<<"first_half: " << std::endl;
    continue_picard = first_half(params,arp,continue_picard);
    //update the effective storativity to use during the next iteration of the first half:
    updateEffectiveStorativity(params,arp);
  }




//get the final T for the halfway point
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) != 0.f)
      arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.wtd_T(x,y), arp.fdepth(x,y), arp.ksat(x,y));
  }




  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      arp.transmissivity(x,y) = ocean_T;
      arp.temp_T(x,y) = ocean_T;
      arp.effective_storativity(x,y) = 1.;
      }
//    else{
      arp.scalar_array_x(x,y) = params.x_partial/arp.effective_storativity(x,y);
      arp.scalar_array_y(x,y) /= arp.effective_storativity(x,y);
 //   }
  }


for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      params.total_loss_to_ocean += arp.wtd(x,y)*arp.cell_area[y];
      arp.wtd(x,y) = 0.;
    }
  }


  //Do the second half of the midpoint method:
  std::cout<<"second_half: " << std::endl;
  second_half(params,arp);


for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      params.total_loss_to_ocean += arp.wtd(x,y)*arp.cell_area[y];
      arp.wtd(x,y) = 0.;
    }
  }

}


void update(Parameters &params, ArrayPack &arp){
  std::cout<<"entering the transient_groundwater module"<<std::endl;
  UpdateCPU(params, arp);
  std::cout<<"leaving the transient_groundwater module"<<std::endl;
}

}
