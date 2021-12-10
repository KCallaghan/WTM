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

      if(continue_picard < 20){
        arp.my_last_wtd(x,y) = (arp.wtd_T(x,y) + arp.wtd_T_iteration(x,y))/2.;
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.my_last_wtd(x,y), arp.fdepth(x,y), arp.ksat(x,y));
        arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
      }
      else{
        if((arp.wtd_T(x,y)>=arp.my_last_wtd(x,y) && arp.wtd_T_iteration(x,y)>=arp.my_last_wtd(x,y)) ||(arp.wtd_T(x,y)<=arp.my_last_wtd(x,y) && arp.wtd_T_iteration(x,y)<=arp.my_last_wtd(x,y))){
          arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);

          if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 2000){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_big;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 200){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_medium;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 10){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_small;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 5){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_tiny;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 1){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_itsy;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 0.3){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_bitsy;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 0.1){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_between;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 0.005){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_spider;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -2000){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) - params.s_big;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -200){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) - params.s_medium;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -10){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_small;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -5){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_tiny;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -1){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_itsy;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -0.3){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_bitsy;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -0.1){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_between;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -0.005){
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_spider;
          }
          else{
            arp.my_last_wtd(x,y) = arp.wtd_T(x,y);
          }
        }
        else{
          double temp = arp.my_last_wtd(x,y);
          arp.my_last_wtd(x,y) = (arp.my_last_wtd(x,y) + arp.my_prev_wtd(x,y))/2.;
          arp.my_prev_wtd(x,y) = temp;
        }

        double new_T = depthIntegratedTransmissivity(arp.my_last_wtd(x,y), arp.fdepth(x,y), arp.ksat(x,y));

        if(arp.nope(x,y) == 1){
	        arp.transmissivity(x,y) = arp.transmissivity(x,y)*0.6 + new_T*0.4;
	      }
      }
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
      b(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);
      guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);
    }
  }

  double entry;
  int main_loc = 0;

  //HALFWAY SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain, provided that they are not ocean cells.
    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
      coefficients.push_back(T(main_loc,main_loc, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(- arp.scalar_array_y_half(x,y)) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(- arp.scalar_array_x_half(x,y)) +1;
      coefficients.push_back(T(main_loc,main_loc, entry));

      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      if(y != params.ncells_y-1 && y!= 0){  // && arp.land_mask(x,y+1) != 0.f){
        //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
        entry = arp.scalar_array_y_half(x,y)*((arp.transmissivity(x,y+1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+1,entry ));
      }
      //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
      if(y != 0 && y != params.ncells_y-1){  // && arp.land_ mask(x,y-1) != 0.f){
        //y may not == 0 since then there is no cell to the west.
        entry = arp.scalar_array_y_half(x,y)*((arp.transmissivity(x,y-1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-1,entry ));
      }
      //Now let's do the North diagonal. Offset by -(ncells_y).
      if(x != 0 && x != params.ncells_x-1){  // && arp.land_mask(x-1,y) != 0.f){
        //x may not equal 0 since then there is no cell to the north.
        entry = arp.scalar_array_x_half(x,y)*((arp.transmissivity(x-1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-params.ncells_y, entry));
      }
      //finally, do the South diagonal, offset by +(ncells_y).
      if(x != params.ncells_x-1 && x!=0){  // && arp.land_mask(x+1,y) != 0.f){
        //we may not be in the final row where there is no cell
        entry = arp.scalar_array_x_half(x,y)*((arp.transmissivity(x+1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+params.ncells_y, entry));
      }
    }
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
  float test_T = 0;

  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){

  //copy result into the wtd_T array:
    if(arp.land_mask(x,y) != 0.f){
      arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
      test_T = depthIntegratedTransmissivity(arp.wtd_T(x,y), arp.fdepth(x,y), arp.ksat(x,y));

      if((1-(std::abs(arp.transmissivity(x,y)/test_T)) < 0.001)  || (std::abs(arp.transmissivity(x,y)-test_T) < 1e-10)){
        arp.nope(x,y) = 0;
      }

      else{
        cell_count += 1;
        arp.nope(x,y) =1;
	    }
    }
  }


  //if(cell_count != 0.){
  //  std::cout<<"The number of cells not equilibrating is "<<cell_count<<". We will continue to picard iteration number "<<picard_number<<std::endl;
  //  picard_number += 1;
  //  return picard_number;
 // }
 // else{
 //   std::cout<<"The number of cells not equilibrating is "<<cell_count<<". Picard iterations terminating."<<std::endl;
  if(picard_number == 10)
    return 0;
  else{
    picard_number += 1;
    return picard_number;
  }
 // }

}






void second_half(const Parameters &params,ArrayPack &arp){

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
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
      guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);
    }
  }

  double entry;
  int main_loc = 0;


  ////SECOND SOLVE
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_loc = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain, provided that they are not ocean cells.
    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
      coefficients_A.push_back(T(main_loc,main_loc, entry));
      coefficients_B.push_back(T(main_loc,main_loc, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
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
    }
  }

  std::cerr<<"finished second set of matrices"<<std::endl;

  //use the triplet vector to populate the matrices, A and B.
  A.setFromTriplets(coefficients_A.begin(),coefficients_A.end());
  B.setFromTriplets(coefficients_B.begin(),coefficients_B.end());


  b = B*b;


 // Apply the second half of the recharge to the water-table depth grid (wtd_T is currently in use). For this, start with the original wtd and add the second half of rech to it.
  //TODO: is this right, or should we start with the halfway-solved wtd_T and add half the rech to it?
 //   #pragma omp parallel for collapse(2)
    for(int y=1;y<params.ncells_y-1;y++)
    for(int x=1;x<params.ncells_x-1; x++){
      if(arp.land_mask(x,y) != 0){          //skip ocean cells
        if(arp.wtd(x,y)>=0){
          b(y+(x*params.ncells_y)) += arp.rech(x,y)/31536000. * params.deltat;
          //if(b(y+(x*params.ncells_y)) <0)
            //b(y+(x*params.ncells_y)) = 0;
        }
        else if (arp.rech(x,y)>0){
          double GW_space = -arp.wtd(x,y) * arp.porosity(x,y);
          if(GW_space > (arp.rech(x,y)/31536000 * params.deltat))
            b(y+(x*params.ncells_y)) += (arp.rech(x,y)/31536000. * params.deltat)/arp.porosity(x,y);
          else{
            double temp = ( (arp.rech(x,y)/31536000. * params.deltat) - GW_space) - arp.wtd(x,y);
            //std::cout<<"x "<<x<<" y "<<y<<" new water table "<<temp<<" old water table "<<
            b(y+(x*params.ncells_y)) += temp;
          }
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
  double x_partial = params.deltat/(params.cellsize_n_s_metres*params.cellsize_n_s_metres);
  float ocean_T = 0.000001 * (1.5 + 25); //some constant for all T values in ocean - TODO look up representative values

std::cout<<"create some needed arrays "<<std::endl;
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      arp.transmissivity(x,y) = ocean_T;
      arp.temp_T(x,y) = ocean_T;
      }
    else{
      arp.scalar_array_x(x,y) = x_partial/arp.porosity(x,y);
      arp.scalar_array_y(x,y) = params.deltat/(arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]*arp.porosity(x,y));
      arp.scalar_array_x_half(x,y) = -arp.scalar_array_x(x,y)/2.;
      arp.scalar_array_y_half(x,y) = -arp.scalar_array_y(x,y)/2.;
    }
  }



 // Apply the first half of the recharge to the water-table depth grid (wtd)
    // Its clone (wtd_T) is used and updated in the Picard iteration
    #pragma omp parallel for collapse(2)
    for(int y=1;y<params.ncells_y-1;y++)
    for(int x=1;x<params.ncells_x-1; x++){
      if(arp.land_mask(x,y) == 0){          //skip ocean cells
        arp.wtd(x,y) = 0;
        //     continue;
      }
      else
        arp.wtd_T(x,y) = add_recharge(params.deltat, arp.rech(x,y)/2., arp.wtd(x,y), arp.land_mask(x,y), arp.porosity(x,y));
        //arp.wtd_T(x,y) = arp.wtd(x,y);
        arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);

    }



  //for (int i=0; i<params.picard_iterations; i++){
   int continue_picard = 1;
   while(continue_picard != 0 ){
    std::cout << "updateTransmissivity: " << std::endl;
    updateTransmissivity(params,arp,continue_picard);
    std::cout<<"first_half: " << std::endl;
    continue_picard = first_half(params,arp,continue_picard);
  }





//get the final T for the halfway point
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) != 0.f)
      arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.wtd_T(x,y), arp.fdepth(x,y), arp.ksat(x,y));
  }

  //Do the second half of the midpoint method:
  std::cout<<"second_half: " << std::endl;
  second_half(params,arp);
}


void update(Parameters &params, ArrayPack &arp){
  std::cout<<"entering the transient_groundwater module"<<std::endl;
  UpdateCPU(params, arp);
  std::cout<<"leaving the transient_groundwater module"<<std::endl;
}

}
