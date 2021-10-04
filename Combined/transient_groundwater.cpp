//#define EIGEN_DONT_PARALLELIZE

#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.
#include <eigen3/Eigen/Core>
//#include <eigen3/Eigen/SPQRSupport>
#include "doctest.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"
#include "mat_mult.cpp"

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

const double FP_ERROR = 1e-4;

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

#define c2d(arr,x,y) arr[(y)*(width)+(x)] //Convert 2d coordinates into 1d coordinates assuming width

struct double2 {
  double x = 0;
  double y = 0;
};


typedef double* d2d_pointer;
typedef float*  f2d_pointer;
typedef double* d1d_pointer;
typedef float*  f1d_pointer;
typedef uint8_t* ui82d_pointer;


struct FanDarcyPack {
  d1d_pointer   cell_area;
  d1d_pointer   cellsize_e_w_metres;
  double        cellsize_n_s_metres;
  d2d_pointer   fdepth;
 // ui82d_pointer land_mask;
  f2d_pointer   land_mask;
  //f2d_pointer   ice_mask;
 
  f2d_pointer   porosity;
  f2d_pointer   topo;
  d2d_pointer   transmissivity;
  d2d_pointer   wtd;
  d2d_pointer   wtd_T;
  d2d_pointer   wtd_changed;
  f2d_pointer   ksat;
  int           width;
};


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
//    std::cout<<"else if 1"<<std::endl;
    return std::max(0.0,fdepth * ksat  * std::exp((wtd_T + shallow)/fdepth));
  } else if(wtd_T > 0){
//    std::cout<<"else if 2 ksat "<<std::endl;
    // If wtd_T is greater than 0, max out rate of groundwater movement
    // as though wtd_T were 0. The surface water will get to move in
    // FillSpillMerge.
    return std::max(0.0,ksat * (0 + shallow + fdepth));
  } else { //Equation S4 from the Fan paper
//    std::cout<<"else"<<std::endl;
    return std::max(0.0,ksat * (wtd_T + shallow + fdepth));  //max because you can't have a negative transmissivity.
  }
}


//update the entire transmissivity array
void updateTransmissivity(
  const Parameters &params,const FanDarcyPack &fdp,ArrayPack &arp
){
  const auto width = fdp.width;
  #pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y-1;y++)
  for(int x=0;x<params.ncells_x-1; x++){
      float ocean_T = 0.000001 * (1.5 + 2.5);  //some constant for all T values in ocean - TODO look up representative values
      if(arp.land_mask(x,y) == 0.f)
        arp.transmissivity(x,y) = ocean_T;
      else
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(c2d(fdp.wtd_T,x,y), c2d(fdp.fdepth,x,y), c2d(fdp.ksat,x,y));
  //    std::cout<<"x "<<x<<" y "<<y<<" transmissivity "<<arp.transmissivity(x,y)<<" fdepth "<<arp.fdepth(x,y)<<" ksat "<<arp.ksat(x,y)<<std::endl;
  }
}





void populateArrays(const Parameters &params,const FanDarcyPack &fdp,ArrayPack &arp){

  std::vector<T> coefficients;
  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);

  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);
  double entry;

    std::cout<<"populate b"<<std::endl;

    std::cout<<"populate triplet"<<std::endl;

    std::cout<<"cells x "<<params.ncells_x<<" cells y "<<params.ncells_y<<std::endl;


  int main_row = 0;
  int main_col = 0;
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    double scalar_portion_x = -params.deltat/(arp.porosity(x,y)*fdp.cellsize_n_s_metres*fdp.cellsize_n_s_metres);
    double scalar_portion_y = -params.deltat/(arp.porosity(x,y)*fdp.cellsize_e_w_metres[y]*fdp.cellsize_e_w_metres[y]);
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
      entry = 1.;//                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.
      coefficients.push_back(T(main_col,main_row, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
      entry = (-2 * arp.transmissivity(x,y))*(scalar_portion_x+scalar_portion_y) +1;//   (-params.deltat/(arp.porosity(x,y)*fdp.cellsize_e_w_metres[y]))*(-4 * arp.transmissivity(x,y)) + 1;
      coefficients.push_back(T(main_col,main_row, entry));

      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
 //       if(!(arp.land_mask(x,y+1) == 0.f)){ //make sure it is not an ocean cell
          if(y==0){
            //y == 0 means we are in the very first column of the domain. There is no cell to the west of the current cell.
            entry = scalar_portion_x*(arp.transmissivity(x,y) + ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));  //ocean_T is used for the T at the non-existant western cell.
            coefficients.push_back(T(main_row,main_col+1,entry ));
          }
          else if(!(y == params.ncells_y-1)){
            //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
            entry = scalar_portion_x*(arp.transmissivity(x,y) + ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));
            coefficients.push_back(T(main_row,main_col+1,entry ));
       //           std::cout<<"east cell x"<<x<<" y "<<y<<" entry "<<entry<<std::endl;

          }
 //       }

        //Next is the West diagonal. Opposite of the East. Located at (i,j-1).
   //     if(!(arp.land_mask(x,y-1) == 0.f)){ //make sure it is not an ocean cell
          if(y== params.ncells_y-1){
            //y is in the final column, there is no cell to the East.
            entry = scalar_portion_x*(arp.transmissivity(x,y) - ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));
            coefficients.push_back(T(main_row,main_col-1,entry ));
          }
          else if(y != 0){
            //y may not == 0 since then there is no cell to the west.
            entry = scalar_portion_x*(arp.transmissivity(x,y) - ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));
            coefficients.push_back(T(main_row,main_col-1,entry ));
          }
  //      }


      //Now let's do the North diagonal. Offset by -(ncells_y).
   //   if(!(arp.land_mask(x-1,y) == 0.f)){ //make sure it is not an ocean cell
        if(x==params.ncells_x-1){
          //we are in the final row of the domain, there is no cell to the South.
          entry = scalar_portion_y*(arp.transmissivity(x,y) - ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
          coefficients.push_back(T(main_row,main_col-params.ncells_y, entry));
        }
        else if(x != 0 ){
          //x may not equal 0 since then there is no cell to the north.
          //also check that it is not an ocean cell.
          entry = scalar_portion_y*(arp.transmissivity(x,y) - ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
          coefficients.push_back(T(main_row,main_col-params.ncells_y, entry));
        }
  //    }

      //finally, do the South diagonal, offset by +(ncells_y).
  //    if(!(arp.land_mask(x+1,y) == 0.f)){ //make sure it is not an ocean cell
        if(x==0){
          //There is no cell to the North of the main cell.
          entry = scalar_portion_y*(arp.transmissivity(x,y) + ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
          coefficients.push_back(T(main_row,main_col+params.ncells_y, entry));
        }
        else if(!(x == params.ncells_x-1)){
          //we may not be in the final row where there is no cell
          entry = scalar_portion_y*(arp.transmissivity(x,y) + ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
          coefficients.push_back(T(main_row,main_col+params.ncells_y, entry));
        }
      }
  //  }
  }

      std::cerr<<"set A"<<std::endl;
//std::cout<<"begin "<<coefficients.begin()<<" end "<<coefficients.end()<<std::endl;
  //use the triplet vector to populate the matrix, A.
  A.setFromTriplets(coefficients.begin(),coefficients.end());
//  std::cout<<"print the matrix"<<std::endl;
//std::cout<< MatrixXd(A)<<std::endl;

//  Eigen::SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
//  vec_x = chol.solve(b);
    std::cerr<<"set solver"<<std::endl;

/*
// Sparse LU solver
Eigen::SparseLU<SpMat, COLAMDOrdering<int> >   solver;
      std::cerr<<"compute"<<std::endl;
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
//solver.analyzePattern(A);
// Compute the numerical factorization
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

//vec_x = lscg.solve(b);
//if(lscg.info()!=Eigen::Success)
 // std::cout<<"solving failed"<<std::endl;
  //vec_x = solver.solve(b);


    std::cerr<<"all donnneeeee"<<std::endl;

//SolverClassName<SparseMatrix<double> > solver;
//solver.compute(A);
//if(solver.info()!=Success) {
//  // decomposition failed
//  std::cout<<"first fail"<<std::endl;
//  return;
//}
//vec_x = solver.solve(b);
//if(solver.info()!=Success) {
//  // solving failed
//  std::cout<<"second fail"<<std::endl;
//  return;
//}


//copy result into the wtd_T array:
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    //if(arp.land_mask(x,y) != 0.f)  //if they are ocean cells
    arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
    //else  //TODO: why does this seem to be necessary? Why am I getting non-zero values in the ocean, ever?
    //  arp.wtd_T(x,y) = 0.;
  }
}




//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void UpdateCPU(const Parameters &params, ArrayPack &arp){
  FanDarcyPack fdp;
  fdp.cell_area           = arp.cell_area.data();
  fdp.cellsize_e_w_metres = arp.cellsize_e_w_metres.data();
  fdp.cellsize_n_s_metres = params.cellsize_n_s_metres;
  fdp.fdepth              = arp.fdepth.data();
  fdp.land_mask           = arp.land_mask.data();
  //fdp.ice_mask            = arp.ice_mask.data();
  fdp.porosity            = arp.porosity.data();
  fdp.topo                = arp.topo.data();
  fdp.transmissivity      = arp.transmissivity.data();
  fdp.width               = arp.fdepth.width();
  fdp.wtd                 = arp.wtd.data(); // External wtd
  fdp.wtd_T               = arp.wtd_T.data(); // Internal wtd for Picard iter
  fdp.wtd_changed         = arp.wtd_changed.data();
  fdp.ksat                = arp.ksat.data();

  Eigen::initParallel();

  // Picard iteration through solver
  // For now, just iterate three times
  int niter = 1;
  for (int i=0; i<niter; i++){
    std::cout << "updateTransmissivity: Iteration " << i+1 << "/" << niter << std::endl;
    updateTransmissivity(params,fdp,arp);
    std::cout<<"populateArrays: Iteration " << i+1 << "/" << niter << std::endl;
    populateArrays(params,fdp,arp);
  }
  // Following these iterations, copy the result into the WTD array
  // >>>> Improve code in future to send results directly to WTD on the
  //      final Picard iteration <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<TODO!
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    //if(arp.land_mask(x,y) != 0.f)  //if they are ocean cells
    arp.wtd(x,y) = arp.wtd_T(x,y);
  }
}


void update(const Parameters &params, ArrayPack &arp){
  std::cout<<"entering the transient_groundwater module"<<std::endl;
  UpdateCPU(params, arp);
  std::cout<<"leaving the transient_groundwater module"<<std::endl;
}


//TEST_CASE("depthIntegratedTransmissivity"){
//  CHECK(depthIntegratedTransmissivity(5  ,100  ,0.01      ) ==doctest::Approx(1.015));
//  CHECK(depthIntegratedTransmissivity(0  ,100  ,0.01      ) ==doctest::Approx(1.015));
//  CHECK(depthIntegratedTransmissivity(-1 ,100  ,0.01      ) ==doctest::Approx(1.005));
//  CHECK(depthIntegratedTransmissivity(-2 ,100  ,0.01      ) ==doctest::Approx(0.995012479192682));
//  CHECK(depthIntegratedTransmissivity(-3 ,100  ,0.1       ) ==doctest::Approx(9.85111939603063));
//  CHECK(depthIntegratedTransmissivity(-4 ,100  ,0.0005    ) ==doctest::Approx(0.048765495601417));
//  CHECK(depthIntegratedTransmissivity(-5 ,100  ,0.0001    ) ==doctest::Approx(0.009656054162576));
//  CHECK(depthIntegratedTransmissivity(-6 ,200  ,0.0001    ) ==doctest::Approx(0.019555024743867));
//  CHECK(depthIntegratedTransmissivity(-7 ,500  ,0.0001    ) ==doctest::Approx(0.049453013938769));
//  CHECK(depthIntegratedTransmissivity(-8 ,1000 ,0.0001    ) ==doctest::Approx(0.099352107930345));
//  CHECK(depthIntegratedTransmissivity(-9 ,10   ,0.0001    ) ==doctest::Approx(0.000472366552741));
//  CHECK(depthIntegratedTransmissivity(-10,500  ,0.0001    ) ==doctest::Approx(0.049157184231746));
//  CHECK(depthIntegratedTransmissivity(-10,10   ,0.0001    ) ==doctest::Approx(0.000427414931949));
//  CHECK(depthIntegratedTransmissivity(-10,1000 ,0.0001    ) ==doctest::Approx(0.099153602286297));
//  CHECK(depthIntegratedTransmissivity(-10,100  ,0.0001    ) ==doctest::Approx(0.009185122844015));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,0.0001    ) ==doctest::Approx(0.029161928740837));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,0.1       ) ==doctest::Approx(29.1619287408366));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,0.5       ) ==doctest::Approx(145.809643704183));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,1         ) ==doctest::Approx(291.619287408366));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,0.000001  ) ==doctest::Approx(0.000291619287408));
//  CHECK(depthIntegratedTransmissivity(-10,300  ,0.0000001 ) ==doctest::Approx(2.91619287408366E-05));
//}
//
//
//
//TEST_CASE("computeNewWTD"){
//  CHECK(computeNewWTD(6000,   5, 0.4, 1000 ) ==doctest::Approx(11));
//  CHECK(computeNewWTD(6000,   1, 0.4, 1000 ) ==doctest::Approx(7));
//  CHECK(computeNewWTD(6000,   1, 0.4, 10000) ==doctest::Approx(1.6));
//  CHECK(computeNewWTD(6000,   1, 0.4, 100  ) ==doctest::Approx(61));
//  CHECK(computeNewWTD(6000,   1, 0.2, 1000 ) ==doctest::Approx(7));
//  CHECK(computeNewWTD(6000,   1, 0.5, 1000 ) ==doctest::Approx(7));
//  CHECK(computeNewWTD(6000,   1, 0.8, 1000 ) ==doctest::Approx(7));
//  CHECK(computeNewWTD(100000, 1, 0.4, 1000 ) ==doctest::Approx(101));
//  CHECK(computeNewWTD(1000,   1, 0.4, 1000 ) ==doctest::Approx(2));
//  CHECK(computeNewWTD(  10,   1, 0.4, 1000 ) ==doctest::Approx(1.01));
//
//  CHECK(computeNewWTD(-6000,   5, 0.4, 1000 ) ==doctest::Approx(-2.5));
//  CHECK(computeNewWTD(-6000,   1, 0.4, 1000 ) ==doctest::Approx(-12.5));
//  CHECK(computeNewWTD(-6000,   1, 0.4, 10000) ==doctest::Approx(0.4));
//  CHECK(computeNewWTD(-6000,   1, 0.4, 100  ) ==doctest::Approx(-147.5));
//  CHECK(computeNewWTD(-6000,   1, 0.2, 1000 ) ==doctest::Approx(-25));
//  CHECK(computeNewWTD(-6000,   1, 0.5, 1000 ) ==doctest::Approx(-10));
//  CHECK(computeNewWTD(-6000,   1, 0.8, 1000 ) ==doctest::Approx(-6.25));
//  CHECK(computeNewWTD(-100000, 1, 0.4, 1000 ) ==doctest::Approx(-247.5));
//  CHECK(computeNewWTD(-1000,   1, 0.4, 1000 ) ==doctest::Approx(0));
//  CHECK(computeNewWTD(  -10,   1, 0.4, 1000 ) ==doctest::Approx(0.99));
//
//  CHECK(computeNewWTD(6000,   -5, 0.4, 1000 ) ==doctest::Approx(4));
//  CHECK(computeNewWTD(6000,   -1, 0.4, 1000 ) ==doctest::Approx(5.6));
//  CHECK(computeNewWTD(6000,   -1, 0.4, 10000) ==doctest::Approx(0.2));
//  CHECK(computeNewWTD(6000,   -1, 0.4, 100  ) ==doctest::Approx(59.6));
//  CHECK(computeNewWTD(6000,   -1, 0.2, 1000 ) ==doctest::Approx(5.8));
//  CHECK(computeNewWTD(6000,   -1, 0.5, 1000 ) ==doctest::Approx(5.5));
//  CHECK(computeNewWTD(6000,   -1, 0.8, 1000 ) ==doctest::Approx(5.2));
//  CHECK(computeNewWTD(100000, -1, 0.4, 1000 ) ==doctest::Approx(99.6));
//  CHECK(computeNewWTD(1000,   -1, 0.4, 1000 ) ==doctest::Approx(0.6));
//  CHECK(computeNewWTD(  10,   -1, 0.4, 1000 ) ==doctest::Approx(-0.975));
//
//  CHECK(computeNewWTD(-6000,   -5, 0.4, 1000 ) ==doctest::Approx(-20));
//  CHECK(computeNewWTD(-6000,   -1, 0.4, 1000 ) ==doctest::Approx(-16));
//  CHECK(computeNewWTD(-6000,   -1, 0.4, 10000) ==doctest::Approx(-2.5));
//  CHECK(computeNewWTD(-6000,   -1, 0.4, 100  ) ==doctest::Approx(-151));
//  CHECK(computeNewWTD(-6000,   -1, 0.2, 1000 ) ==doctest::Approx(-31));
//  CHECK(computeNewWTD(-6000,   -1, 0.5, 1000 ) ==doctest::Approx(-13));
//  CHECK(computeNewWTD(-6000,   -1, 0.8, 1000 ) ==doctest::Approx(-8.5));
//  CHECK(computeNewWTD(-100000, -1, 0.4, 1000 ) ==doctest::Approx(-251));
//  CHECK(computeNewWTD(-1000,   -1, 0.4, 1000 ) ==doctest::Approx(-3.5));
//  CHECK(computeNewWTD(  -10,   -1, 0.4, 1000 ) ==doctest::Approx(-1.025));
//
//}


//TODO: Example array test case
//TEST_CASE("calculateWaterVolume"){
//  f2d porosity = {{1,1,1},{2,2,2},{3,3,3}};
//  f2d fdepth = {{1,1,1},{2,2,2},{3,3,3}};
//  std::vector<double> cell_area = {1,2,3};
//
//  FanDarcyPack fdp;
//  fdp.porosity  = porosity.data();
//  fdp.fdepth    = porosity.data();
//  fdp.cell_area = cell_area.data();
//  fdp.width     = porosity.width();
//
//  //NOTE: Use SUBCASE if you need to reconstruct the arrays multiple times (e.g. if they get changed)
//  //See: https://github.com/onqtam/doctest/blob/master/doc/markdown/tutorial.md#test-cases-and-subcases
//  calculateWaterVolume(3, 4, 5, 1, 1, 2, 2, fdp);
//}

}
