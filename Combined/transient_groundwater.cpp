#include "doctest.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"
#include "mat_mult.cpp"

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
  f2d_pointer   ice_mask;
 
  f2d_pointer   porosity;
  f2d_pointer   topo;
  d2d_pointer   transmissivity;
  d2d_pointer   wtd;
  d2d_pointer   wtd_changed;
  f2d_pointer   ksat;
  int           width;
};


double depthIntegratedTransmissivity(
  const double wtd,
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
  } else if(wtd < -shallow){ // Equation S6 from the Fan paper
    return std::max(0.0,fdepth * ksat  * std::exp((wtd + shallow)/fdepth));
  } else if(wtd > 0){
    // If wtd is greater than 0, max out rate of groundwater movement
    // as though wtd were 0. The surface water will get to move in
    // FillSpillMerge.
    return std::max(0.0,ksat * (0 + shallow + fdepth));
  } else { //Equation S4 from the Fan paper
    return std::max(0.0,ksat * (wtd + shallow + fdepth));  //max because you can't have a negative transmissivity.
  }
}


//update the entire transmissivity array
void updateTransmissivity(
  const Parameters &params,const FanDarcyPack &fdp,ArrayPack &arp
){
  const auto width = fdp.width;
  #pragma omp parallel for collapse(2)
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)          //skip ocean cells
      continue;
    arp.transmissivity(x,y) = depthIntegratedTransmissivity(c2d(fdp.wtd,x,y), c2d(fdp.fdepth,x,y), c2d(fdp.ksat,x,y));
  }
}



void populateArrays(const Parameters &params,const FanDarcyPack &fdp,ArrayPack &arp){
  for(int y=0;y<params.ncells_y;y++)
  for(int x=0;x<params.ncells_x; x++){  //not sure why it is ncells-1. Should it just be ncells?

    //populate the diagonal matrix
    arp.diagonal_matrix(x+(y*params.ncells_y),x+(y*params.ncells_y)) = (-params.deltat/(arp.porosity(x,y)*fdp.cellsize_e_w_metres[y]))*(-4 * arp.transmissivity(x,y)) + 1;  //which direction of cellsize am I supposed to use here?
    if(x!=0 && x != params.ncells_x-1){ //what happens in the corner cells? I can't add data because we need both T_i+1 and T_i-1??
      arp.diagonal_matrix(x+(y*params.ncells_y)-1,x+(y*params.ncells_y)) = (-params.deltat/(arp.porosity(x-1,y)*fdp.cellsize_e_w_metres[y]))*(arp.transmissivity(x,y) - ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
      arp.diagonal_matrix(x+(y*params.ncells_y)+1,x+(y*params.ncells_y)) = (-params.deltat/(arp.porosity(x+1,y)*fdp.cellsize_e_w_metres[y]))*(arp.transmissivity(x,y) + ((arp.transmissivity(x+1,y)-arp.transmissivity(x-1,y))/4));
    }
    if(x!=0 && x!=1 && y!=0 &&y!=params.ncells_y-1) //what happens in the corner cells? I can't add data because we need both T_i+1 and T_i-1??
      arp.diagonal_matrix(x+(y*params.ncells_y)-2,x+(y*params.ncells_y)) = (-params.deltat/(arp.porosity(x,y-1)*fdp.cellsize_n_s_metres))*(arp.transmissivity(x,y) - ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));
    if(x != params.ncells_x-2 && x!=params.ncells_x-1 && y!=0 &&y!=params.ncells_y-1)
      arp.diagonal_matrix(x+(y*params.ncells_y)+2,x+(y*params.ncells_y)) = (-params.deltat/(arp.porosity(x,y+1)*fdp.cellsize_n_s_metres))*(arp.transmissivity(x,y) + ((arp.transmissivity(x,y+1)-arp.transmissivity(x,y-1))/4));


    //populate the starting guess vector for the wtd
    arp.wtd_1D(x+(y*params.ncells_y)) = arp.wtd(x,y);
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
  fdp.ice_mask            = arp.ice_mask.data();
  fdp.porosity            = arp.porosity.data();
  fdp.topo                = arp.topo.data();
  fdp.transmissivity      = arp.transmissivity.data();
  fdp.width               = arp.fdepth.width();
  fdp.wtd                 = arp.wtd.data();
  fdp.wtd_changed         = arp.wtd_changed.data();
  fdp.ksat                = arp.ksat.data();

std::cout<<"updateTransmissivity"<<std::endl;
updateTransmissivity(params,fdp,arp);

std::cout<<"populateArrays"<<std::endl;
populateArrays(params,fdp,arp);

std::cout<<"mat_mult"<<std::endl;
mat_mult(params.ncells_x*params.ncells_y,params.ncells_x*params.ncells_y, params.ncells_x*params.ncells_y,1, arp.diagonal_matrix,arp.wtd_1D,arp.wtd_1D_out);


std::cout<<"update wtd"<<std::endl;
  // Once all the changes are known, update the WTD everywhere with the
  // difference array
  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=0; y<params.ncells_y; y++)
  for(int32_t x=0; x<params.ncells_x; x++){
    // Update the whole wtd array at once.
    // This is the new water table after groundwater has moved
    // for delta_t seconds.
    arp.wtd(x,y) = arp.wtd_1D_out(x+y*params.ncells_y);;
  }
arp.wtd.printAll("answer so far - after");
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
