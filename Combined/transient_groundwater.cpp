#include "add_recharge.hpp"
#include "transient_groundwater.hpp"
#include "update_effective_storativity.hpp"

using namespace Eigen;

typedef Eigen::SparseMatrix<double,RowMajor> SpMat; // declares a row-major sparse matrix type of double
typedef Eigen::Triplet<double> T;  // used to populate the matrices

constexpr double solver_tolerance_value = 0.00001;
constexpr double seconds_in_a_year = 31536000.;

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
  Parameters &params,ArrayPack &arp
){
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) != 0.f){
        double my_last_wtd = (arp.wtd_T(x,y) + arp.wtd_T_iteration(x,y))/2.;
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(my_last_wtd, arp.fdepth(x,y), static_cast<double>(arp.ksat(x,y)));
    }
  }
}



// The midpoint method consists of two main steps.
// In the first step, we compute water tables at the time that is *half* of the full time-step.
// To do so, we use an implicit backward-difference Euler method.
// Using these half-time water tables, we compute the new transmissivity values
// That will be used in the second step of the midpoint method below.
void first_half(const Parameters &params,ArrayPack &arp){

  std::vector<T> coefficients;
  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x*params.ncells_y);
  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);


  //We need to solve the vector-matrix equation Ax=b.
  //b consists of the current head values (i.e. water table depth + topography)
  //We also populate a 'guess', which consists of a water table (wtd_T) that in later iterations
  //has already been modified for changing transmissivity closer to the final answer.
  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);
    b(y+(x*params.ncells_y)) = arp.wtd(x,y) + static_cast<double>(arp.topo(x,y));  //wtd is 0 in ocean cells and topo is 0 in ocean cells, so no need to differentiate between ocean vs land.
    guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + static_cast<double>(arp.topo(x,y));
  }

  double entry;
  int main_loc = 0; //the cell to populate in the A matrix.

  //HALFWAY SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=1;x<params.ncells_x-1; x++)
  for(int y=1;y<params.ncells_y-1; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain.
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(arp.scalar_array_y(x,y)/(arp.effective_storativity(x,y))) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(params.x_partial/(2*arp.effective_storativity(x,y))) +1;
      coefficients.push_back(T(main_loc,main_loc, entry));

      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      entry = -arp.scalar_array_y(x,y)/(arp.effective_storativity(x,y))*((arp.transmissivity(x,y+1) + arp.transmissivity(x,y))/2.);
      coefficients.push_back(T(main_loc,main_loc+1,entry ));

      //Next is the West diagonal. Opposite of the East. Located at (i,j-1).
      entry = -arp.scalar_array_y(x,y)/(arp.effective_storativity(x,y))*((arp.transmissivity(x,y-1) + arp.transmissivity(x,y))/2.);
      coefficients.push_back(T(main_loc,main_loc-1,entry ));

      //Now let's do the North diagonal. Offset by -(ncells_y).
      entry = -params.x_partial/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x-1,y) + arp.transmissivity(x,y))/2.);
      coefficients.push_back(T(main_loc,main_loc-params.ncells_y, entry));

      //finally, do the South diagonal, offset by +(ncells_y).
      entry = -params.x_partial/(2*arp.effective_storativity(x,y))*((arp.transmissivity(x+1,y) + arp.transmissivity(x,y))/2.);
      coefficients.push_back(T(main_loc,main_loc+params.ncells_y, entry));
  }

  //use the triplet vector to populate the matrix, A.
  A.setFromTriplets(coefficients.begin(),coefficients.end());


  // Biconjugate gradient solver with guess
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(solver_tolerance_value);
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
    arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - static_cast<double>(arp.topo(x,y));
  }
}



// In the second step of the midpoint method, we use the transmissivity values
// obtained after the calculation in the first step above.
// We then compute the new water table depths after a full time-step has passed.
void second_half(Parameters &params,ArrayPack &arp){

  std::vector<T> coefficients_A;
  std::vector<T> coefficients_B;
  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);
  SpMat B(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);

  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x*params.ncells_y);

  //We need to solve the vector-matrix equation Ax=Bb.
  //b consists of the current head values (i.e. water table depth + topography)
  //We also populate a 'guess', which consists of a water table (wtd_T) that
  //has already been modified for changing transmissivity closer to the final answer.
  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    b(y+(x*params.ncells_y)) = arp.original_wtd(x,y) + static_cast<double>(arp.topo(x,y));  //original wtd is 0 in ocean cells and topo is 0 in ocean cells, so no need to differentiate between ocean vs land.
    guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + static_cast<double>(arp.topo(x,y));
  }

  double entry;
  int main_loc = 0;

  ////SECOND SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=1;x<params.ncells_x-1; x++)
  for(int y=1;y<params.ncells_y-1; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain.
      entry =  1. - (arp.scalar_array_x(x,y)*(arp.transmissivity(x-1,y)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x+1,y)/4.)) - (arp.scalar_array_y(x,y)*(arp.transmissivity(x,y-1)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x,y+1)/4.));
      coefficients_B.push_back(T(main_loc,main_loc, entry));
      entry =  1. + (arp.scalar_array_x(x,y)*(arp.transmissivity(x-1,y)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x+1,y)/4.)) + (arp.scalar_array_y(x,y)*(arp.transmissivity(x,y-1)/4. + arp.transmissivity(x,y)/2. + arp.transmissivity(x,y+1)/4.));
      coefficients_A.push_back(T(main_loc,main_loc, entry));

     //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      entry = arp.scalar_array_y(x,y)*(arp.transmissivity(x,y) + arp.transmissivity(x,y+1))/4.;
      coefficients_B.push_back(T(main_loc,main_loc+1,entry ));
      coefficients_A.push_back(T(main_loc,main_loc+1,-entry ));

      //Next is the West diagonal. Opposite of the East. Located at (i,j-1).
      entry = arp.scalar_array_y(x,y)*(arp.transmissivity(x,y) + arp.transmissivity(x,y-1))/4.;
      coefficients_B.push_back(T(main_loc,main_loc-1,entry ));
      coefficients_A.push_back(T(main_loc,main_loc-1,-entry ));

      //Now let's do the North diagonal. Offset by -(ncells_y).
      entry = arp.scalar_array_x(x,y)*(arp.transmissivity(x,y) + arp.transmissivity(x-1,y))/4.;
      coefficients_B.push_back(T(main_loc,main_loc-params.ncells_y,entry ));
      coefficients_A.push_back(T(main_loc,main_loc-params.ncells_y,-entry ));

      //finally, do the South diagonal, offset by +(ncells_y).
      entry = arp.scalar_array_x(x,y)*(arp.transmissivity(x,y) + arp.transmissivity(x+1,y))/4.;
      coefficients_B.push_back(T(main_loc,main_loc+params.ncells_y,entry ));
      coefficients_A.push_back(T(main_loc,main_loc+params.ncells_y,-entry ));
  }

  std::cerr<<"finished second set of matrices"<<std::endl;

  //use the triplet vector to populate the matrices, A and B.
  A.setFromTriplets(coefficients_A.begin(),coefficients_A.end());
  B.setFromTriplets(coefficients_B.begin(),coefficients_B.end());


  b = B*b;
  guess = B*guess;


   //Apply the recharge to the water-table depth grid. The full amount of recharge is added to b, which was created based on original_wtd.
   //Recharge is added here because of the form that the equation takes - can't add it prior to doing the B*b multiplication.
    for(int y=0;y<params.ncells_y;y++)
    for(int x=0;x<params.ncells_x; x++){
      if(arp.land_mask(x,y) == 1){        //skip ocean cells
        double rech_change = arp.rech(x,y)/seconds_in_a_year * params.deltat;
        params.total_added_recharge += rech_change*arp.cell_area[y];
        if(arp.original_wtd(x,y) >= 0){          //there was surface water, so recharge may be negative
          b(y+(x*params.ncells_y)) += rech_change;
          guess(y+(x*params.ncells_y)) += rech_change;
          if(arp.original_wtd(x,y) + rech_change <0){
            double temp = -(arp.original_wtd(x,y) + rech_change);
            b(y+(x*params.ncells_y)) += temp;
            guess(y+(x*params.ncells_y)) += temp;
            params.total_added_recharge += temp*arp.cell_area[y];  //in this scenario, there has been a negative amount of recharge, so temp here will be positive and remove the extra amount that was spuriously subtracted.
          }
        }
        else if (rech_change>0){        //when there is no surface water, only positive changes in recharge are allowed
          double GW_space = -arp.original_wtd(x,y) * static_cast<double>(arp.porosity(x,y));
          if(GW_space > rech_change){
            b(y+(x*params.ncells_y)) += rech_change/static_cast<double>(arp.porosity(x,y));
            guess(y+(x*params.ncells_y)) += rech_change/static_cast<double>(arp.porosity(x,y));
          }
          else{
            double temp = ( rech_change - GW_space) - arp.original_wtd(x,y);
            b(y+(x*params.ncells_y)) += temp;
            guess(y+(x*params.ncells_y)) += temp;
          }
        }
      }
    }


  // Biconjugate gradient solver with guess
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(solver_tolerance_value);
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
  //copy result into the wtd array:
    arp.wtd(x,y) = vec_x(y+(x*params.ncells_y)) - static_cast<double>(arp.topo(x,y));
  }
  std::cerr<<"finished assigning the new wtd"<<std::endl;
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
  double ocean_T = 0.000001 * (1.5 + 25.); //some constant for all T values in ocean - TODO look up representative values
  //this assumes a ksat of 0.000001 and an e-folding depth of 25. 1.5 is a standard value based on the shallow depths to which soil textures are known.

  std::cout<<"create some needed arrays "<<std::endl;

  //no pragma because we're editing params.total_loss_to_ocean
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      params.total_loss_to_ocean += arp.wtd(x,y)*arp.cell_area[y];
      arp.wtd(x,y) = 0.;
    }
  }

  //Apply the first half of the recharge to the water-table depth grid (wtd)
  //Its clone (wtd_T) is used and updated in the Picard iteration
  //also set the starting porosity
  //set the scalar arrays for x and y directions
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      arp.transmissivity(x,y) = ocean_T;
      arp.wtd_T(x,y) = 0.;
      arp.original_wtd(x,y) = 0.;
      arp.wtd_T_iteration(x,y) = 0.;
      arp.effective_storativity(x,y) = 1.;
    }
    else{
      arp.original_wtd(x,y) = arp.wtd(x,y);
      arp.wtd(x,y) = add_recharge(params.deltat/2., arp.rech(x,y), arp.wtd(x,y), arp.land_mask(x,y), static_cast<double>(arp.porosity(x,y))); //use regular porosity for adding recharge since this checks for underground space within add_recharge.
      arp.wtd_T(x,y) = arp.wtd(x,y);
      arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);

      if(arp.original_wtd(x,y)>0)
        arp.effective_storativity(x,y) = 1;
      else
        arp.effective_storativity(x,y) = static_cast<double>(arp.porosity(x,y));
    }
    arp.scalar_array_y(x,y) = params.deltat/(2.*arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]);
  }


   for(int continue_picard = 0;continue_picard<3;continue_picard++){
    std::cout << "updateTransmissivity: " << std::endl;
    updateTransmissivity(params,arp);
    std::cout<<"first_half: " << std::endl;
    first_half(params,arp);
    //update the effective storativity to use during the next iteration of the first half:
    for(int y=0;y<params.ncells_y;y++)
    for(int x=0;x<params.ncells_x; x++){
      if(arp.land_mask(x,y) != 0.f)
        arp.effective_storativity(x,y) = updateEffectiveStorativity(arp.original_wtd(x,y),arp.wtd_T(x,y), arp.porosity(x,y), arp.effective_storativity(x,y));
    }
  }


  //get the final T for the halfway point
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) != 0.f){
      arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.wtd_T(x,y), arp.fdepth(x,y), static_cast<double>(arp.ksat(x,y)));
    }
    arp.scalar_array_x(x,y) = params.x_partial/arp.effective_storativity(x,y);
    arp.scalar_array_y(x,y) = params.deltat/(arp.effective_storativity(x,y)*arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]);
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
