//#define EIGEN_DONT_PARALLELIZE

#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.
#include <eigen3/Eigen/Core>
#include "doctest.h"
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
  const Parameters &params,ArrayPack &arp
){
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y-1;y++)
  for(int x=0;x<params.ncells_x-1; x++){
      float ocean_T = 0.000001 * (1.5 + 25);  //some constant for all T values in ocean - TODO look up representative values
      if(arp.land_mask(x,y) == 0.f)
        arp.transmissivity(x,y) = ocean_T;
      else
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.wtd_T(x,y), arp.fdepth(x,y), arp.ksat(x,y));
  }
}



int populateArrays(const Parameters &params,ArrayPack &arp, int picard_number){

  std::vector<T> coefficients;
  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);

  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);
  double entry;

    std::cout<<"populate b and triplet"<<std::endl;

    std::cout<<"cells x "<<params.ncells_x<<" cells y "<<params.ncells_y<<std::endl;


  int main_row = 0;
  int main_col = 0;
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    double scalar_portion_x = -params.deltat/(arp.porosity(x,y)*params.cellsize_n_s_metres*params.cellsize_n_s_metres);
    double scalar_portion_y = -params.deltat/(arp.porosity(x,y)*arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]);
    //The row and column that the current cell will be stored in in matrix A.
    //This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
    //All of the N,E,S,W directions should be in the same row, but the column will differ.
    main_row = y+(x*params.ncells_y);
    main_col = y+(x*params.ncells_y);
    //start with the central diagonal, which contains the info for the current cell:
    //This diagonal will be populated for all cells in the domain, provided that they are not ocean cells.
    if(arp.land_mask(x,y) == 0.f){  //if they are ocean cells
      //populate the known vector b. This is the current wtd, which is the 'guess' that we are using to get our answer, x.
      b(y+(x*params.ncells_y)) = 0.;
      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
      coefficients.push_back(T(main_row,main_col, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
      if(y== 608 && x ==1595)
        std::cout<<"x "<<x<<" y "<<y<<" wtd "<<arp.wtd(x,y)<<" b "<<b(y+(x*params.ncells_y))<<" wtd_T "<<arp.wtd_T(x,y) <<std::endl;
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(- scalar_portion_x) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(- scalar_portion_y) +1;
      coefficients.push_back(T(main_row,main_col, entry));

      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
          if(!(y == params.ncells_y-1)){  // && arp.land_mask(x,y+1) != 0.f){
            //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
            entry = scalar_portion_x*((3./8.)*arp.transmissivity(x,y+1) + arp.transmissivity(x,y)/4. + (3./8.)*arp.transmissivity(x,y-1));
            coefficients.push_back(T(main_row,main_col+1,entry ));
          }

        //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
          if(y != 0){  // && arp.land_ mask(x,y-1) != 0.f){
            //y may not == 0 since then there is no cell to the west.
            entry = scalar_portion_x*((1./8.)*arp.transmissivity(x,y+1) + (3./4.)*arp.transmissivity(x,y) + (1./8.)*arp.transmissivity(x,y-1));
            coefficients.push_back(T(main_row,main_col-1,entry ));
          }


      //Now let's do the North diagonal. Offset by -(ncells_y).
        if(x != 0 ){  // && arp.land_mask(x-1,y) != 0.f){
          //x may not equal 0 since then there is no cell to the north.
          entry = scalar_portion_y*((1./8.)*arp.transmissivity(x+1,y) + (3./4.)*arp.transmissivity(x,y) + (1./8.)*arp.transmissivity(x-1,y));
          coefficients.push_back(T(main_row,main_col-params.ncells_y, entry));
        }

      //finally, do the South diagonal, offset by +(ncells_y).
        if(!(x == params.ncells_x-1)){  // && arp.land_mask(x+1,y) != 0.f){
           //we may not be in the final row where there is no cell
          entry = scalar_portion_y*((3./8.)*arp.transmissivity(x+1,y) + arp.transmissivity(x,y)/4. + (3./8.)*arp.transmissivity(x-1,y));
          coefficients.push_back(T(main_row,main_col+params.ncells_y, entry));
        }
      }
  }

      std::cerr<<"set A"<<std::endl;

  //use the triplet vector to populate the matrix, A.
  A.setFromTriplets(coefficients.begin(),coefficients.end());

    std::cerr<<"set solver"<<std::endl;

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
//NOTE: we cannot use the Eigen:IncompleteLUT preconditioner, because its implementation is serial. Using it means that BiCGSTAB will not run in parallel. It is faster without.
std::cout<<"compute"<<std::endl;
solver.compute(A);
//solver.setTolerance(1e-15);
assert(solver.info()==Eigen::Success);
std::cout<<"solve"<<std::endl;
vec_x = solver.solveWithGuess(b, b);  // guess = b;

std::cout<<"set the new wtd_T values "<<std::endl;
//copy result into the wtd_T array:
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
  }




int cell_count = 0;
int total_cells = params.ncells_x * params.ncells_y;
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    if(std::abs(arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)) > 0.1)
        cell_count += 1;
    arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);
  }


if(cell_count < total_cells/100.){
  std::cout<<"The number of cells not equilibrating is "<<cell_count<<", which is less than 1 percent of the total cells. Picard iterations terminating."<<std::endl;
  return 0;
}
else{
  picard_number += 1;
  std::cout<<"The number of cells not equilibrating is "<<cell_count<<", out of a total of "<<total_cells<<" cells. We will continue to picard iteration number "<<picard_number<<std::endl;
  return picard_number;
}



}



//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void UpdateCPU(const Parameters &params, ArrayPack &arp){

  Eigen::initParallel();
  omp_set_num_threads(8);
  Eigen::setNbThreads(8);

  // Picard iteration through solver

  //for (int i=0; i<params.picard_iterations; i++){
   int continue_picard = 1;
   while(continue_picard != 0 ){
    std::cout << "updateTransmissivity: " << std::endl;
    updateTransmissivity(params,arp);
    std::cout<<"populateArrays: " << std::endl;
    continue_picard = populateArrays(params,arp,continue_picard);
  }
  // Following these iterations, copy  the result into the WTD array
  // >>>> Improve code in future to send results directly to WTD on the
  //      final Picard iteration <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TODO!
  std::cout<<"update the wtd to the new set of values "<<std::endl;
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
      arp.wtd(x,y) = arp.wtd_T(x,y);
  }
}


void update(const Parameters &params, ArrayPack &arp){
  std::cout<<"entering the transient_groundwater module"<<std::endl;
  UpdateCPU(params, arp);
  std::cout<<"leaving the transient_groundwater module"<<std::endl;
}

}