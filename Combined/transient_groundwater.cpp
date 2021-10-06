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

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
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



void populateArrays(const Parameters &params,ArrayPack &arp){

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
      coefficients.push_back(T(main_col,main_row, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
      if(y== 608 && x ==1595)
        std::cout<<"x "<<x<<" y "<<y<<" wtd "<<arp.wtd(x,y)<<" b "<<b(y+(x*params.ncells_y))<<" wtd_T "<<arp.wtd_T(x,y) <<std::endl;
      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(- scalar_portion_x) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(- scalar_portion_y) +1;
      coefficients.push_back(T(main_col,main_row, entry));

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

/*
// Biconjugate gradient solver
//Eigen::BiCGSTAB<SpMat> solver;
Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
      std::cerr<<"compute"<<std::endl;
solver.analyzePattern(A);
solver.compute(A);
      std::cerr<<"check"<<std::endl;
    std::cerr<<"solve"<<std::endl;
vec_x = solver.solve(b);
*/


// Biconjugate gradient solver with guess
// Set up the guess -- same as last time's levels (b)
// or topography (arp.topo)
// Just guessing it is b now; commenting out!
/*
for(int x=0;x<params.ncells_x; x++)
for(int y=0;y<params.ncells_y; y++){
    //vec_x(y+(x*params.ncells_y)) = b(y+(x*params.ncells_y));
    //vec_x(y+(x*params.ncells_y)) = arp.topo(x,y);
    vec_x(y+(x*params.ncells_y)) = 500.;
}
*/
// Solver
//Eigen::BiCGSTAB<SpMat> solver;
Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
      std::cerr<<"compute"<<std::endl;
//solver.analyzePattern(A);
solver.compute(A);
      std::cerr<<"check"<<std::endl;
    std::cerr<<"solve"<<std::endl;
// guess = b; use first line otherwise & uncomment the above
//vec_x = solver.solveWithGuess(b, vec_x);
vec_x = solver.solveWithGuess(b, b);


/*
 // Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int> > solver;
//Eigen::BiCGSTAB<SpMat> solver;
//Eigen::PartialPivLU<SpMat> solver;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
//Eigen::LeastSquaresConjugateGradient<SpMat> solver;
      std::cerr<<"compute"<<std::endl;
solver.compute(A);
      std::cerr<<"check"<<std::endl;
//solver.factorize(A);

assert(solver.info()==Eigen::Success);
//Use the factors to solve the linear system

 //   Eigen::SPQR<SpMat> lscg(A);
//lscg.compute(A);
//if(lscg.info()!=Eigen::Success)
//  std::cerr<<"failed 222 "<<lscg.info()<<std::endl;
//Eigen::SimplicialLDLT
//lscg.setTolerance(1e-14);

 // if(solver.info() != Eigen::Success) {
    // decomposition failed
   // std::cout<<"first failed"<<std::endl;
  //  return;
 // }

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
vec_x = solver.solve(b);
*/

/*
// Biconjugate gradient solver
//Eigen::BiCGSTAB<SpMat> solver;
Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
      std::cerr<<"compute"<<std::endl;
solver.analyzePattern(A);
solver.compute(A);
      std::cerr<<"check"<<std::endl;
    std::cerr<<"solve"<<std::endl;
vec_x = solver.solve(b);
*/


// Biconjugate gradient solver with guess
// Set up the guess -- same as last time's levels (b)
// or topography (arp.topo)
// Just guessing it is b now; commenting out!
/*
for(int x=0;x<params.ncells_x; x++)
for(int y=0;y<params.ncells_y; y++){
    //vec_x(y+(x*params.ncells_y)) = b(y+(x*params.ncells_y));
    //vec_x(y+(x*params.ncells_y)) = arp.topo(x,y);
    vec_x(y+(x*params.ncells_y)) = 500.;
}
*/
// Solver
//Eigen::BiCGSTAB<SpMat> solver;
Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
solver.compute(A);
assert(solver.info()==Eigen::Success);
vec_x = solver.solveWithGuess(b, b);  // guess = b;

std::cout<<"set the new wtd_T values "<<std::endl;
//copy result into the wtd_T array:
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
  }
}




//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void UpdateCPU(const Parameters &params, ArrayPack &arp){

   Eigen::initParallel();

  // Picard iteration through solver

  for (int i=0; i<params.picard_iterations; i++){
    std::cout << "updateTransmissivity: Iteration " << i+1 << "/" << params.picard_iterations << std::endl;
    updateTransmissivity(params,arp);
    std::cout<<"populateArrays: Iteration " << i+1 << "/" << params.picard_iterations << std::endl;
    populateArrays(params,arp);
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
