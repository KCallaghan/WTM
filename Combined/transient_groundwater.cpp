//#define EIGEN_DONT_PARALLELIZE

#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.
#include <eigen3/Eigen/Core>
#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"
#include "mat_mult.cpp"

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
	  if(continue_picard>100)
		std::cout<<"x "<<x<<" y "<<y<<" old T "<<arp.transmissivity(x,y)<<" new T "<<new_T<<std::endl;
          arp.transmissivity(x,y) = arp.transmissivity(x,y)*0.6 + new_T*0.4;
	  if(continue_picard>100)
	      std::cout<<"x "<<x<<" y "<<y<<" updated T is "<<arp.transmissivity(x,y)<<std::endl;
	}
      }
    }
  }
}



int populateArrays(const Parameters &params,ArrayPack &arp, int picard_number){

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
      //populate the known vector b. This is the current wtd, which is the 'guess' that we are using to get our answer, x.
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
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(- arp.scalar_array_y(x,y)) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(- arp.scalar_array_x(x,y)) +1;

      coefficients.push_back(T(main_loc,main_loc, entry));


      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
      if(y != params.ncells_y-1 && y!= 0){  // && arp.land_mask(x,y+1) != 0.f){
        //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
        entry = arp.scalar_array_y(x,y)*((arp.transmissivity(x,y+1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+1,entry ));
      }
      //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
      if(y != 0 && y != params.ncells_y-1){  // && arp.land_ mask(x,y-1) != 0.f){
        //y may not == 0 since then there is no cell to the west.
        entry = arp.scalar_array_y(x,y)*((arp.transmissivity(x,y-1) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-1,entry ));
      }
      //Now let's do the North diagonal. Offset by -(ncells_y).
      if(x != 0 && x != params.ncells_x-1){  // && arp.land_mask(x-1,y) != 0.f){
        //x may not equal 0 since then there is no cell to the north.
        entry = arp.scalar_array_x(x,y)*((arp.transmissivity(x-1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc-params.ncells_y, entry));
      }
      //finally, do the South diagonal, offset by +(ncells_y).
      if(x != params.ncells_x-1 && x!=0){  // && arp.land_mask(x+1,y) != 0.f){
        //we may not be in the final row where there is no cell
        entry = arp.scalar_array_x(x,y)*((arp.transmissivity(x+1,y) + arp.transmissivity(x,y))/2.);
        coefficients.push_back(T(main_loc,main_loc+params.ncells_y, entry));
      }
    }
  }


  //use the triplet vector to populate the matrix, A.
  A.setFromTriplets(coefficients.begin(),coefficients.end());


/*
// Sparse LU solver
Eigen::SparseLU<SpMat, COLAMDOrdering <int> >   solver;
      std::cerr<<"compute"<<std:: endl;
solver.compute(A);
      std::cerr<<"check"<<std::endl;
    std::cerr<<"solve"<<std::endl;
vec_x = solver.solve(b);
*/

/*
// UMFPACK -- not installed / set up
Eigen::UmfPackLU<SpMat>   solver;
      std::cerr<<"compute"<<std::endl;
solver.compute(A);
      std::cerr<<"check"<<std::endl;
    std::cerr<<"solve"<<std::endl;
vec_x = solver.solve (b);
*/


//#pragma omp parallel
//    {  printf("Hello World from thread = %d\n", omp_get_thread_num()); }


//int nthreads = Eigen::nbThreads( );
//std::cout << "THREADS = " << nthreads <<std::ends;


// Biconjugate gradient solver with guess
// Set up the guess -- same as last time's levels (b)
Eigen::BiCGSTAB<SpMat> solver;//, Eigen::IncompleteLUT<double> > solver;
solver.setTolerance(0.00001);
//Eigen::LeastSquaresConjugateGradient<SpMat> solver;
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

int ten_m = 0;
int one_m = 0;
int ten_cm = 0;
int one_cm = 0;
int one_mm = 0;
int tenth_mm = 0;
int less = 0;
int cells_to_adjust = 0;

  //#pragma omp parallel for collapse(2)
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
        cells_to_adjust +=1;
	  if(picard_number>100)
		  std::cout<<"x "<<x<<" y "<<y<<" wtd_T "<<arp.wtd_T(x,y)<<" wtd_T_iteration "<<arp.wtd_T_iteration(x,y)<<" transmissivity "<<arp.transmissivity(x,y)<<" test_T "<<test_T<<" fdepth "<<arp.fdepth(x,y)<<" last wtd "<<arp.my_last_wtd(x,y)<<std::endl;
      }





//Print stats on how close we are:
if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 10) && arp.land_mask(x,y) !=0.f)
  ten_m +=1;
else if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 1)&& arp.land_mask(x,y) !=0.f)
  one_m += 1;
else if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 0.1)&& arp.land_mask(x,y) !=0.f)
  ten_cm +=1;
else if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 0.01)&& arp.land_mask(x,y) != 0.f){
  one_cm +=1;
//  if (picard_number>100)
  //  std::cout<<"x "<<x<<" y "<<y<<std::endl;
}
else if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 0.001)&& arp.land_mask(x,y) != 0.f)
  one_mm +=1;
else if((std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 0.0001)&& arp.land_mask(x,y) != 0.f)
  tenth_mm +=1;
else if (arp.land_mask(x,y) != 0.f)
  less +=1;



    }
  }


std::cout<<"here are the stats on how close we are: "<<std::endl;
std::cout<<ten_m<<" cells are more than 10 m off; "<<std::endl;
std::cout<<one_m<<" cells are more than 1 m off; "<<std::endl;
std::cout<<ten_cm<<" cells are more than 10 cm off; "<<std::endl;
std::cout<<one_cm<<" cells are more than 1 cm off; "<<std::endl;
std::cout<<one_mm<<" cells are more than 1 mm off; "<<std::endl;
std::cout<<tenth_mm<<" cells are more than 1/10 mm off; "<<std::endl;
std::cout<<less<<" cells are closer than 1/10 mm; "<<std::endl;
std::cout<<"cells_to_adjust: "<<cells_to_adjust<<std::endl;




//  if(picard_number == 200){
//    std::cout<<"The number of cells not equilibrating is "<<cell_count<<". 100 picard iterations, so we will terminate."<<std::endl;
//    return 0;
//  }
  if(cell_count != 0.){
    std::cout<<"The number of cells not equilibrating is "<<cell_count<<". We will continue to picard iteration number "<<picard_number<<std::endl;
    picard_number += 1;
    return picard_number;
  }
  else{
    std::cout<<"The number of cells not equilibrating is "<<cell_count<<". Picard iterations terminating."<<std::endl;
    return 0;
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
  double x_partial = -params.deltat/(params.cellsize_n_s_metres*params.cellsize_n_s_metres);
  float ocean_T = 0.000001 * (1.5 + 25); //some constant for all T values in ocean - TODO look up representative values

  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    if(arp.land_mask(x,y) == 0.f){
      arp.transmissivity(x,y) = ocean_T;
      }
    else{
      arp.scalar_array_x(x,y) = x_partial/arp.porosity(x,y);
      arp.scalar_array_y(x,y) = -params.deltat/(arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]*arp.porosity(x,y));
    }
  }



  //for (int i=0; i<params.picard_iterations; i++){
   int continue_picard = 1;
   while(continue_picard != 0 ){
    std::cout << "updateTransmissivity: " << std::endl;
    updateTransmissivity(params,arp,continue_picard);
    std::cout<<"populateArrays: " << std::endl;
    continue_picard = populateArrays(params,arp,continue_picard);
  }
  // Following these iterations, copy  the result into the WTD array
  // >>>> Improve code in future to send results directly to WTD on the
  //      final Picard iteration <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TODO!
  std::cout<<"update the wtd to the new set of values "<<std::endl;

  #pragma omp parallel for collapse(2)
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
      arp.wtd(x,y) = arp.wtd_T(x,y);
  }
}


void update(Parameters &params, ArrayPack &arp){
  std::cout<<"entering the transient_groundwater module"<<std::endl;
  UpdateCPU(params, arp);
  std::cout<<"leaving the transient_groundwater module"<<std::endl;
}

}
