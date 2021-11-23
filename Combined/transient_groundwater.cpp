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
    return std::max(0.0,fdepth * ksat  * std::exp((wtd_T + shallow)/(fdepth*5)));
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
int big_count = 0;
int medium_count = 0;
int small_count = 0;
int tiny_count = 0;
int itsy_count = 0;
int bitsy_count = 0;
int tracking_count = 0;
int fluctuating_count = 0;
int ocean_count = 0;
int total_count = 0;
int spider_count = 0;
 //     if(continue_picard % 100 == 0){
 //       params.s_big = params.s_big/3.;
 //       params.s_medium = params.s_medium/3.;
 //       params.s_small = params.s_small/3.;
 //       params.s_tiny = params.s_tiny/3.;
 //       params.s_itsy =params.s_itsy/3.;
 //       params.s_bitsy = params.s_bitsy/3.;
 //       std::cout<<"the new values are "<<params.s_big<<" med "<<params.s_medium<<" small "<<params.s_small<<" tiny "<<params.s_tiny<<std::endl;
 //     }
  //#pragma omp parallel for collapse(2)
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){
    total_count +=1;

      float ocean_T = 0.000001 * (1.5 + 25);  //some constant for all T values in ocean - TODO look up representative values


      if(arp.land_mask(x,y) == 0.f){
        if(x==88&&y==70)
          std::cout<<"it's an ocean cell"<<std::endl;
        arp.transmissivity(x,y) = ocean_T;
        arp.my_last_wtd(x,y) = arp.wtd_T(x,y);
        ocean_count +=1;
      }
      else if(continue_picard < 20){
        arp.my_last_wtd(x,y) = (arp.wtd_T(x,y) + arp.wtd_T_iteration(x,y))/2.;
        arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.my_last_wtd(x,y), arp.fdepth(x,y), arp.ksat(x,y));
        arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);

      }
      else{
       // if(arp.wtd_T(x,y)<=0 && arp.my_last_wtd(x,y) > 0){
       //   arp.my_last_wtd(x,y) = 0;
       // }
        if((arp.wtd_T(x,y)>=arp.my_last_wtd(x,y) && arp.wtd_T_iteration(x,y)>=arp.my_last_wtd(x,y)) ||(arp.wtd_T(x,y)<=arp.my_last_wtd(x,y) && arp.wtd_T_iteration(x,y)<=arp.my_last_wtd(x,y))){

          if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 2000){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_big;
            big_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 200){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_medium;
            medium_count += 1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 10){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) + params.s_small;
            small_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 5){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_tiny;
            tiny_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 1){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_itsy;
            itsy_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 0.3){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_bitsy;
            bitsy_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) > 0.005){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) +params.s_spider;
            spider_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -2000){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) - params.s_big;
            big_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -200){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) - params.s_medium;
            medium_count += 1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -10){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_small;
            small_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -5){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_tiny;
            tiny_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -1){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_itsy;
            itsy_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -0.3){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_bitsy;
            bitsy_count +=1;
          }
          else if(arp.wtd_T(x,y)- arp.my_last_wtd(x,y) < -0.005){
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.my_last_wtd(x,y) -params.s_spider;
            spider_count +=1;
          }
          else{
            arp.my_prev_wtd(x,y) = arp.my_last_wtd(x,y);
            arp.my_last_wtd(x,y) = arp.wtd_T(x,y);
            tracking_count +=1;
          }
        }
        else{
          double temp = arp.my_last_wtd(x,y);
          arp.my_last_wtd(x,y) = (arp.my_last_wtd(x,y) + arp.my_prev_wtd(x,y))/2.;
          arp.my_prev_wtd(x,y) = temp;
          fluctuating_count +=1;

        }


//if(arp.wtd_T(x,y)<=0 && arp.my_last_wtd(x,y) > 0)
 // arp.my_last_wtd(x,y) = arp.wtd_T(x,y)/10.;


        arp.transmissivity(x,y) = depthIntegratedTransmissivity(arp.my_last_wtd(x,y), arp.fdepth(x,y), arp.ksat(x,y));
      }


  }


std::cout<<"big_count "<<big_count<<" medium_count "<<medium_count<<" small_count "<<small_count<<" tiny_count "<<tiny_count<<" itsy_count "<<itsy_count<<" bitsy_count "<<bitsy_count<<" spider_count "<<spider_count<<" tracking_count "<<tracking_count<<" fluctuating_count "<<fluctuating_count<<" ocean_count "<<ocean_count<<" total_count "<<total_count<<std::endl;

}



int populateArrays(const Parameters &params,ArrayPack &arp, int picard_number){

  std::vector<T> coefficients;
  Eigen::VectorXd b(params.ncells_x*params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x*params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x*params.ncells_y);

  SpMat A(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y);
  double entry;

    std::cout<<"populate b and triplet"<<std::endl;

    std::cout<<"cells x "<<params.ncells_x<<" cells y "<<params.ncells_y<<std::endl;

  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    arp.wtd_T_iteration(x,y) = arp.wtd_T(x,y);
  }


  int main_row = 0;
  int main_col = 0;
  //populate the coefficients triplet vector. This should have row index, column index, value of what is needed in the final matrix A.
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y; y++){
    double scalar_portion_x = -params.deltat/(arp.porosity(x,y)*params.cellsize_n_s_metres*params.cellsize_n_s_metres);
    double scalar_portion_y = -params.deltat/(arp.porosity(x,y)*arp.cellsize_e_w_metres[y]*arp.cellsize_e_w_metres[y]);

//TODO: Do I have my usage of scalar_portion_x and scalar_portion_y the right way around? Getting confused if I want the distance travelled or the cross-sectional area that the T is felt across.

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
            guess(y+(x*params.ncells_y)) = 0.;

      entry = 1.;                   //then they should not change (wtd should be 0 both before and after) so only the centre diagonal is populated, and it is = 1.

      coefficients.push_back(T(main_row,main_col, entry));
    }
    else{  //land cells, so we have an actual value here and we should consider the neighbouring cells.
      b(y+(x*params.ncells_y)) = arp.wtd(x,y) + arp.topo(x,y);
            guess(y+(x*params.ncells_y)) = arp.wtd_T(x,y) + arp.topo(x,y);

      entry =  (arp.transmissivity(x,y-1)/2 + arp.transmissivity(x,y) + arp.transmissivity(x,y+1)/2)*(- scalar_portion_y) + (arp.transmissivity(x-1,y)/2 + arp.transmissivity(x,y) + arp.transmissivity(x+1,y)/2)*(- scalar_portion_x) +1;
      //entry = arp.transmissivity(x,y)*-2*(-scalar_portion_y - scalar_portion_x ) + 1;

      coefficients.push_back(T(main_row,main_col, entry));


      //Now do the East diagonal. Because C++ is row-major, the East location is at (i,j+1).
          if(y != params.ncells_y-1 && y!= 0){  // && arp.land_mask(x,y+1) != 0.f){
            //y may not == params.ncells_y-1, since this would be the eastern-most cell in the domain. There is no neighbour to the east.
            //entry = scalar_portion_y*((3./8.)*arp.transmissivity(x,y+1) + arp.transmissivity(x,y)/4. + (3./8.)*arp.transmissivity(x,y-1));
       //     entry = scalar_portion_y*((3./8.)*arp.transmissivity(x,y+1) + arp.transmissivity(x,y)/2. + (1./8.)*arp.transmissivity(x,y-1));
            entry = scalar_portion_y*((arp.transmissivity(x,y+1) + arp.transmissivity(x,y))/2.);
              //entry = scalar_portion_y*(arp.transmissivity(x,y) + (arp.transmissivity(x,y+1) + arp.transmissivity(x,y-1))/4.);

            coefficients.push_back(T(main_row,main_col+1,entry ));
          }

        //Next is the West diagonal. Opposite of the East. Lo cated at (i,j-1).
          if(y != 0 && y != params.ncells_y-1){  // && arp.land_ mask(x,y-1) != 0.f){
            //y may not == 0 since then there is no cell to the west.
           // entry = scalar_portion_y*((1./8.)*arp.transmissivity(x,y+1) + (1./2.)*arp.transmissivity(x,y) + (3./8.)*arp.transmissivity(x,y-1));
            entry = scalar_portion_y*((arp.transmissivity(x,y-1) + arp.transmissivity(x,y))/2.);
            //entry = scalar_portion_y*(arp.transmissivity(x,y) - (arp.transmissivity(x,y+1) + arp.transmissivity(x,y-1))/4.);
            coefficients.push_back(T(main_row,main_col-1,entry ));
          }


      //Now let's do the North diagonal. Offset by -(ncells_y).
        if(x != 0 && x != params.ncells_x-1){  // && arp.land_mask(x-1,y) != 0.f){
          //x may not equal 0 since then there is no cell to the north.
          //entry = scalar_portion_x*((1./8.)*arp.transmissivity(x+1,y) + (1./2.)*arp.transmissivity(x,y) + (3./8.)*arp.transmissivity(x-1,y));
          entry = scalar_portion_x*((arp.transmissivity(x-1,y) + arp.transmissivity(x,y))/2.);
          //entry = scalar_portion_x*(arp.transmissivity(x,y) - (arp.transmissivity(x+1,y) + arp.transmissivity(x-1,y))/4.);
          coefficients.push_back(T(main_row,main_col-params.ncells_y, entry));
        }

      //finally, do the South diagonal, offset by +(ncells_y).
        if(x != params.ncells_x-1 && x!=0){  // && arp.land_mask(x+1,y) != 0.f){
           //we may not be in the final row where there is no cell
          //entry = scalar_portion_x*((3./8.)*arp.transmissivity(x+1,y) + arp.transmissivity(x,y)/2. + (1./8.)*arp.transmissivity(x-1,y));
          entry = scalar_portion_x*((arp.transmissivity(x+1,y) + arp.transmissivity(x,y))/2.);
          //entry = scalar_portion_x*(arp.transmissivity(x,y) + (arp.transmissivity(x+1,y) + arp.transmissivity(x-1,y))/4.);
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
solver.setTolerance(0.000001);
//Eigen::LeastSquaresConjugateGradient<SpMat> solver;
//NOTE: we cannot use the Eigen:IncompleteLUT preconditioner, because its implementation is serial. Using it means that BiCGSTAB will not run in parallel. It is faster without.
std::cout<<"compute"<<std::endl;
solver.compute(A);

if(solver.info()!=Success) {
  std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ decomposition failed!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
}

//solver.setTolerance(1e-15);
assert(solver.info()==Eigen::Success);
std::cout<<"solve"<<std::endl;
vec_x = solver.solveWithGuess(b, guess);  // guess = b;

if(solver.info()!=Success) {
  std::cout<<"########################################################################### solving failed!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
}


  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

std::cout<<"set the new wtd_T values "<<std::endl;
//copy result into the wtd_T array:
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    arp.wtd_T(x,y) = vec_x(y+(x*params.ncells_y)) - arp.topo(x,y);
  }







int cell_count = 0;
int total_cells = params.ncells_x * params.ncells_y;
int land_cells = 0;
int land_eq = 0;
  for(int x=0;x<params.ncells_x; x++)
  for(int y=0;y<params.ncells_y;y++){
    if(arp.land_mask(x,y) != 0.f){
      land_cells +=1;
    }


 //         if(x==68&&y==34)
  //          if((x==89||x==88||x==90)&&(y==71||y==70||y==72) )
  //    std::cout<<"x "<<x<<" y "<<y<<" wtd_T "<<arp.wtd_T(x,y)<<" wtd_T_iteration "<<arp.wtd_T_iteration(x,y)<<" my wtd "<<arp.my_last_wtd(x,y)<<" prev was "<<arp.my_prev_wtd(x,y)<<" transmissivity "<<arp.transmissivity(x,y)<<std::endl;

    if((std::abs(arp.my_last_wtd(x,y) - arp.wtd_T(x,y)) > 0.001)){
      cell_count += 1;
      arp.nope(x,y)=1;
          if(arp.land_mask(x,y) != 0.f){
            land_eq +=1;
          }

 //     if(picard_number == 10)
   //     std::cout<<"x "<<x<<" y "<<y<<" diff "<<arp.wtd_T_iteration(x,y) - arp.wtd_T(x,y)<<std::endl;
    }
    else{
      arp.nope(x,y)=0;
    }
  }


if(cell_count < land_cells/100.){
  std::cout<<"The number of cells not equilibrating is "<<cell_count<<", which is less than 1 percent of the total cells. Picard iterations terminating."<<std::endl;
  std::cout<<"There are a total of "<<land_cells<<" land cells and of these, "<<land_eq<<" are not equilibrating."<<std::endl;
  return 0;
}
else{
  picard_number += 1;
  std::cout<<"The number of cells not equilibrating is "<<cell_count<<", out of a total of "<<total_cells<<" cells. We will continue to picard iteration number "<<picard_number<<std::endl;
  std::cout<<"There are a total of "<<land_cells<<" land cells and of these, "<<land_eq<<" are not equilibrating."<<std::endl;
  return picard_number;
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