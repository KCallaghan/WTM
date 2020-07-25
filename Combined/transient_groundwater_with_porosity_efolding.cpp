#include "doctest.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"

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
  f2d_pointer   fdepth;
  ui82d_pointer land_mask;
  f2d_pointer   porosity;
  f2d_pointer   topo;
  f2d_pointer   transmissivity;
  f2d_pointer   wtd;
  f2d_pointer   wtd_changed;
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
    return fdepth * ksat * std::exp((wtd + shallow)/fdepth);
  } else if(wtd > 0){
    // If wtd is greater than 0, max out rate of groundwater movement
    // as though wtd were 0. The surface water will get to move in
    // FillSpillMerge.
    return ksat * (0 + shallow + fdepth);
  } else { //Equation S4 from the Fan paper
    return ksat * (wtd + shallow + fdepth);
  }
}



void updateTransmissivity(
  const int x,
  const int y,
  std::array<double,5> &local_wtd,
  const FanDarcyPack &fdp
){
  const auto width = fdp.width;

  c2d(fdp.transmissivity,x,y) = depthIntegratedTransmissivity(
      local_wtd[0], c2d(fdp.fdepth,x,y), c2d(fdp.ksat,x,y));
  c2d(fdp.transmissivity,x,y+1) = depthIntegratedTransmissivity(
      local_wtd[1], c2d(fdp.fdepth,x,y+1), c2d(fdp.ksat,x,y+1));
  c2d(fdp.transmissivity,x,y-1) = depthIntegratedTransmissivity(
      local_wtd[2], c2d(fdp.fdepth,x,y-1), c2d(fdp.ksat,x,y-1));
  c2d(fdp.transmissivity,x+1,y) = depthIntegratedTransmissivity(
      local_wtd[3], c2d(fdp.fdepth,x+1,y), c2d(fdp.ksat,x+1,y));
  c2d(fdp.transmissivity,x-1,y) = depthIntegratedTransmissivity(
      local_wtd[4], c2d(fdp.fdepth,x-1,y), c2d(fdp.ksat,x-1,y));
}


double vonNeumannStability(
  const double Dmax,

  const double cell_width,
  const double cell_height
){
  constexpr double SAFETY = 2;
  return
      (                 std::pow(cell_width, 2) * std::pow(cell_height, 2)) /
      (SAFETY * Dmax * (std::pow(cell_width, 2) + std::pow(cell_height, 2)));
}


/**
 * @brief Returns the maximum stable time step with a 2x factor of safety
 * @details Uses a 2D diffusion von Neumann stability analysis using the
 *          "worst case" scenario highest transmissivity, smallest cellsize, 
 *          combined with a porosity-based amplification factor.
 */
double computeMaxStableTimeStep(
  const int x,
  const int y,
  const FanDarcyPack &fdp
){
  const int width = fdp.width;

  // Transmissivity is an effective diffusivity
  // Use the highest transmissivity for this time-step calculation

  const std::array<double,5> Tarray = {
    c2d(fdp.transmissivity, x  , y),
    c2d(fdp.transmissivity, x  , y+1),
    c2d(fdp.transmissivity, x  , y-1),
    c2d(fdp.transmissivity, x-1, y  ),
    c2d(fdp.transmissivity, x+1, y  )
  };


  const std::array<double,3> sizeArray = {
    fdp.cellsize_e_w_metres[y],
    fdp.cellsize_e_w_metres[y+1],
    fdp.cellsize_e_w_metres[y-1]
  };

  const auto Dmax = *std::max_element(Tarray.begin(), Tarray.end());

  const auto sizeMin = *std::min_element(sizeArray.begin(), sizeArray.end());
  const auto dt_max_diffusion_basic = vonNeumannStability(Dmax, sizeMin, fdp.cellsize_n_s_metres);


const double e_folding_depths_center = c2d(fdp.wtd, x  , y  )/c2d(fdp.fdepth, x  , y  );
const double e_folding_depths_N = c2d(fdp.wtd, x  , y+1  )/c2d(fdp.fdepth, x  , y+1  );
const double e_folding_depths_S = c2d(fdp.wtd, x  , y-1  )/c2d(fdp.fdepth, x  , y-1  );
const double e_folding_depths_E = c2d(fdp.wtd, x+1  , y  )/c2d(fdp.fdepth, x+1  , y  );
const double e_folding_depths_W = c2d(fdp.wtd, x-1  , y  )/c2d(fdp.fdepth, x-1  , y  );

//TODO: should we calculate this using the lowest head of all 5 cells, since this is where head could potentially lower to?
//TODO: surface porosities should be set to min 0.2. 
 double porosityCenter = c2d(fdp.porosity, x  , y  )/exp(e_folding_depths_center);
  double porosityN = c2d(fdp.porosity, x  , y+1  )/exp(e_folding_depths_N);
  double porosityS = c2d(fdp.porosity, x  , y-1  )/exp(e_folding_depths_S);
  double porosityE = c2d(fdp.porosity, x+1  , y  )/exp(e_folding_depths_E);
  double porosityW = c2d(fdp.porosity, x-1  , y  )/exp(e_folding_depths_W);

 if(c2d(fdp.wtd, x  , y  ) > -1.5)
    porosityCenter = c2d(fdp.porosity, x  , y  );
  if(c2d(fdp.wtd, x  , y+1  ) > -1.5)
    porosityN = c2d(fdp.porosity, x  , y+1  );
  if(c2d(fdp.wtd, x  , y-1  ) > -1.5)
    porosityS = c2d(fdp.porosity, x  , y-1  );
  if(c2d(fdp.wtd, x+1  , y  ) > -1.5)
    porosityE = c2d(fdp.porosity, x+1  , y  );
  if(c2d(fdp.wtd, x-1  , y  ) > -1.5)
    porosityW = c2d(fdp.porosity, x-1  , y  );

  // Now let's add in the porosity differences
  // Porosity differences will AMPLIFY the WTD changes.
  // Let us choose a time step that is also based on the LOWEST porosity
  // (highest WTD change per water volume transfer)
  auto PhiMin =             porosityCenter;
  PhiMin = std::min(PhiMin, porosityN);
  PhiMin = std::min(PhiMin, porosityS);
  PhiMin = std::min(PhiMin, porosityE);
  PhiMin = std::min(PhiMin, porosityW);

  // Porosity is a linear amplifier of WTD change, and it amplifies change
  // in both the giving and receiving cells.
  // Amplification goes as 1 / phi.
  // We need to consider this going up and down, so the amplification becomes
  // 2/Phi, and therefore we multiply our stable dt by Phi/2
  const auto dt_max_diffusion_withPorosity = dt_max_diffusion_basic * PhiMin / 2.;
  // In order to avoid operating at the very maximum time step possible,
  // we apply a factor of safety of 2

  return dt_max_diffusion_withPorosity/2.;
}


double computeNewWTD(
  double volume_change,
  const double initial_wtd,
  const double fdepth,
  const double porosity_at_surface,
  const double cell_area
){
  // DEFINE VARIABLES
  const double FP_ERROR = 0.00001;
  // Total Groundwater Storage Capacity
  const auto total_gw_storage_capacity = cell_area * porosity_at_surface * fdepth;
  // New water-table depth: to return
  double final_wtd;

  // IF THE CELL STARTS WITH SURFACE WATER
  if ( initial_wtd >= 0 ){
    // Check if final water-table depth will also be > 0
    // If it is, then either the cell is gaining water, 
    // or it is losing water, but it is losing less than the 
    // amount of surface water it has. 
    const auto initial_surface_water_volume = initial_wtd * cell_area;
    // Final water volume with respect to the land surface.
    // Positive: above surface
    // Negative: below surface
    auto Vfinal_surface = initial_surface_water_volume + volume_change;
    if ( Vfinal_surface >= 0 ){
      final_wtd = Vfinal_surface/cell_area;  //the cell still has surface water.
    }
    // Otherwise, if WTD <= 0 (i.e. subsurface)
    else{
      //if it's a losing cell, we need to make sure that it's not going to try
      //to lose more water than it has.

      assert(-Vfinal_surface <= total_gw_storage_capacity + FP_ERROR);
  //    if(-Vfinal_surface > total_gw_storage_capacity)
   //     Vfinal_surface = -total_gw_storage_capacity + FP_ERROR; 
        //The volume moving may not be exactly equal to the volume available, or the 
        //logarithm below will be undefined, so we add FP_ERROR to make it just slightly smaller.
      final_wtd = fdepth * std::log1p(Vfinal_surface/total_gw_storage_capacity);
    }
  }
  // OTHERWISE, IF THE CELL STARTS WITH WATER TABLE IN THE SUBSURFACE
  else{
    // Remaining subsurface dry pore volume
    const auto initial_remaining_subsurface_pore_volume = total_gw_storage_capacity * -std::expm1(initial_wtd/fdepth);
    // Final water volume with respect to the land surface.
    // Positive: above surface
    // Negative: below surface
    const auto Vfinal_surface = volume_change - initial_remaining_subsurface_pore_volume;
    // If the volume stays below the surface for the whole time
    if ( Vfinal_surface < 0 ){
      //if it's a losing cell, we need to make sure that it's not going to try
      //to lose more water than it has.
      if(volume_change < 0){
        assert(-volume_change <= (total_gw_storage_capacity - initial_remaining_subsurface_pore_volume + FP_ERROR));
    //    if(-volume_change > (total_gw_storage_capacity - initial_remaining_subsurface_pore_volume))
     //     volume_change = -(total_gw_storage_capacity - initial_remaining_subsurface_pore_volume - FP_ERROR); 
        //The volume moving may not be exactly equal to the volume available, or the 
        //logarithm below will be undefined, so we add FP_ERROR to make it just slightly smaller.
      }
      final_wtd = fdepth * std::log( volume_change/total_gw_storage_capacity
                                     + std::exp(initial_wtd/fdepth) );

    }
    // Otherwise, if the water goes above the surface
    else{
      final_wtd = Vfinal_surface/cell_area;
    }
  }

  return final_wtd;
}


/**
 * @brief Calculates water-table depth change at a cell (and associated
 *        surrounding cells) and updates the class variables associated
 *        with these.
 */
void computeWTDchangeAtCell(
  const int x,
  const int y,
  const double dt,
  std::array<double,5> &local_wtd,
  const FanDarcyPack &fdp
){
  const int width = fdp.width;
  // Update WTD change

  // We do this instead of using a staggered grid to approx. double CPU time
  // in exchange for using less memory.

  // First, compute elevation head at center cell and all neighbours
  // This equals topography plus the water table depth
  // (positive upwards; negative if water table is below Earth surface)
  const double headCenter = c2d(fdp.topo, x  , y  ) + local_wtd[0];
  const double headN      = c2d(fdp.topo, x  , y+1) + local_wtd[1];
  const double headS      = c2d(fdp.topo, x  , y-1) + local_wtd[2];
  const double headE      = c2d(fdp.topo, x+1, y  ) + local_wtd[3];
  const double headW      = c2d(fdp.topo, x-1, y  ) + local_wtd[4];

  const double transmissivityN = std::min(c2d(fdp.transmissivity, x, y),c2d(fdp.transmissivity, x  , y+1));
  const double transmissivityS = std::min(c2d(fdp.transmissivity, x, y),c2d(fdp.transmissivity, x  , y-1));
  const double transmissivityE = std::min(c2d(fdp.transmissivity, x, y),c2d(fdp.transmissivity, x+1  , y));
  const double transmissivityW = std::min(c2d(fdp.transmissivity, x, y),c2d(fdp.transmissivity, x-1  , y));

  // Then, compute the discharges
  // Define Q such that + if center cell gaining, - if center cell losing
  const double QN = transmissivityN * (headN - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QS = transmissivityS * (headS - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QE = transmissivityE * (headE - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;
  const double QW = transmissivityW * (headW - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;

  // Sum discharges and divide by the time step to compute the water volume
  // added to the center cell
  const double dVolume = (QN + QS + QW + QE) * dt;// / (fdp.cellsize_n_s_metres * fdp.cellsize_e_w_metres[y]);


const double wtd_before = local_wtd[0];
const double wtdN_before = local_wtd[1];
const double wtdS_before = local_wtd[2];
const double wtdE_before = local_wtd[3];
const double wtdW_before = local_wtd[4];



  // Update the cell's WTD
  local_wtd[0] = computeNewWTD( dVolume, local_wtd[0], c2d(fdp.fdepth, x, y), c2d(fdp.porosity, x, y), fdp.cell_area[y] );

  // For local dynamic time stepping (consider switching to global later),
  // update the neighboring cell WTDs
  local_wtd[1] = computeNewWTD( -QN*dt, local_wtd[1], c2d(fdp.fdepth, x, y+1), c2d(fdp.porosity, x, y+1), fdp.cell_area[y+1] );
  local_wtd[2] = computeNewWTD( -QS*dt, local_wtd[2], c2d(fdp.fdepth, x, y-1), c2d(fdp.porosity, x, y-1), fdp.cell_area[y-1] );
  local_wtd[3] = computeNewWTD( -QE*dt, local_wtd[3], c2d(fdp.fdepth, x+1, y), c2d(fdp.porosity, x+1, y), fdp.cell_area[y] );
  local_wtd[4] = computeNewWTD( -QW*dt, local_wtd[4], c2d(fdp.fdepth, x-1, y), c2d(fdp.porosity, x-1, y), fdp.cell_area[y] );

}



/**
 * @brief Updates the wtd_depth_total array at a cell(x,y) using the
 * pre-set time step and dynamic time stepping within this as needed.
 */
double updateCell(
  const int x,
  const int y,
  double time_remaining,
  const FanDarcyPack &fdp
){
  const auto width = fdp.width;

  // Skip ocean cells
  if(c2d(fdp.land_mask,x,y) != 1)
    return 0;                     //the coastline represents a boundary condition where water table is at the land surface

  // Runs functions to compute time steps and update WTD for the center cell
  // and its neighbours until the outer time step has been completed

  // Initial water-table depths, prior to updating
  std::array<double, 5> local_wtd = {
    c2d(fdp.wtd, x  , y  ),
    c2d(fdp.wtd, x  , y+1),
    c2d(fdp.wtd, x  , y-1),
    c2d(fdp.wtd, x+1, y  ),
    c2d(fdp.wtd, x-1, y  )
  };

  // Update water-table depths using dynamic time stepping
  double dt_inner;
  while (time_remaining > 0){
    updateTransmissivity(x, y, local_wtd, fdp);

    const double max_stable_time_step = computeMaxStableTimeStep(x, y, fdp);
    // Choose the inner-loop time step
    if(time_remaining <= max_stable_time_step){
      dt_inner = time_remaining;
    }
    else{
      dt_inner = max_stable_time_step;
    }

    computeWTDchangeAtCell(x, y, dt_inner, local_wtd, fdp);
    time_remaining -= dt_inner;
  }

  // When exiting loop, the local_wtd[0] variable holds the final water-table depth
  return local_wtd[0];
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
  fdp.porosity            = arp.porosity.data();
  fdp.topo                = arp.topo.data();
  fdp.transmissivity      = arp.transmissivity.data();
  fdp.width               = arp.fdepth.width();
  fdp.wtd                 = arp.wtd.data();
  fdp.wtd_changed         = arp.wtd_changed.data();
  fdp.ksat                = arp.ksat.data();


  // Updates water-table depth grid over one time step
  #pragma omp parallel for collapse(2) default(none) shared(arp,params,fdp)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    arp.wtd_changed(x,y) = updateCell(x, y, params.deltat, fdp);
  }

  // Once all the changes are known, update the WTD everywhere with the
  // difference array
  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    // Update the whole wtd array at once.
    // This is the new water table after groundwater has moved
    // for delta_t seconds.
    arp.wtd(x,y) = arp.wtd_changed(x,y);;
  }
}



void update(const Parameters &params, ArrayPack &arp){
  UpdateCPU(params, arp);
}


TEST_CASE("depthIntegratedTransmissivity"){
  CHECK(depthIntegratedTransmissivity(5  ,100  ,0.01      ) ==doctest::Approx(1.015));
  CHECK(depthIntegratedTransmissivity(0  ,100  ,0.01      ) ==doctest::Approx(1.015));
  CHECK(depthIntegratedTransmissivity(-1 ,100  ,0.01      ) ==doctest::Approx(1.005));
  CHECK(depthIntegratedTransmissivity(-2 ,100  ,0.01      ) ==doctest::Approx(0.995012479192682));
  CHECK(depthIntegratedTransmissivity(-3 ,100  ,0.1       ) ==doctest::Approx(9.85111939603063));
  CHECK(depthIntegratedTransmissivity(-4 ,100  ,0.0005    ) ==doctest::Approx(0.048765495601417));
  CHECK(depthIntegratedTransmissivity(-5 ,100  ,0.0001    ) ==doctest::Approx(0.009656054162576));
  CHECK(depthIntegratedTransmissivity(-6 ,200  ,0.0001    ) ==doctest::Approx(0.019555024743867));
  CHECK(depthIntegratedTransmissivity(-7 ,500  ,0.0001    ) ==doctest::Approx(0.049453013938769));
  CHECK(depthIntegratedTransmissivity(-8 ,1000 ,0.0001    ) ==doctest::Approx(0.099352107930345));
  CHECK(depthIntegratedTransmissivity(-9 ,10   ,0.0001    ) ==doctest::Approx(0.000472366552741));
  CHECK(depthIntegratedTransmissivity(-10,500  ,0.0001    ) ==doctest::Approx(0.049157184231746));
  CHECK(depthIntegratedTransmissivity(-10,10   ,0.0001    ) ==doctest::Approx(0.000427414931949));
  CHECK(depthIntegratedTransmissivity(-10,1000 ,0.0001    ) ==doctest::Approx(0.099153602286297));
  CHECK(depthIntegratedTransmissivity(-10,100  ,0.0001    ) ==doctest::Approx(0.009185122844015));
  CHECK(depthIntegratedTransmissivity(-10,300  ,0.0001    ) ==doctest::Approx(0.029161928740837));
  CHECK(depthIntegratedTransmissivity(-10,300  ,0.1       ) ==doctest::Approx(29.1619287408366));
  CHECK(depthIntegratedTransmissivity(-10,300  ,0.5       ) ==doctest::Approx(145.809643704183));
  CHECK(depthIntegratedTransmissivity(-10,300  ,1         ) ==doctest::Approx(291.619287408366));
  CHECK(depthIntegratedTransmissivity(-10,300  ,0.000001  ) ==doctest::Approx(0.000291619287408));
  CHECK(depthIntegratedTransmissivity(-10,300  ,0.0000001 ) ==doctest::Approx(2.91619287408366E-05));
}



TEST_CASE("computeNewWTD"){
  CHECK(computeNewWTD(487.660600188348  ,-10  ,300, 0.1, 1000      ) ==doctest::Approx(-5));
  CHECK(computeNewWTD(877.789080339026  ,-10  ,300, 0.2, 900       ) ==doctest::Approx(-5));
  CHECK(computeNewWTD(746.262468812392  ,-3   ,300, 0.5, 500       ) ==doctest::Approx(0));
  CHECK(computeNewWTD(1194.01552781598  ,-2   ,300, 0.8, 1500      ) ==doctest::Approx(-1));
  CHECK(computeNewWTD(215.318057232321  ,-100 ,300, 0.3, 1000      ) ==doctest::Approx(-99));
  CHECK(computeNewWTD(1391.76019394264  ,-10  ,100, 0.3, 1000      ) ==doctest::Approx(-5));
  CHECK(computeNewWTD(1488.79363305426  ,-10  ,1000,0.3, 1000      ) ==doctest::Approx(-5));
  CHECK(computeNewWTD(715.953655623573  ,-10  ,10,  0.3, 1000      ) ==doctest::Approx(-5));

  CHECK(computeNewWTD(4000              ,1    ,300, 0.1, 1000      ) ==doctest::Approx(5));
  CHECK(computeNewWTD(4500              ,0    ,300, 0.2, 900       ) ==doctest::Approx(5));  
  CHECK(computeNewWTD(5000              ,10   ,300, 0.5, 500       ) ==doctest::Approx(20));  
  CHECK(computeNewWTD(15000             ,0    ,300, 0.8, 1500      ) ==doctest::Approx(10));  

  CHECK(computeNewWTD(5983.51698553982  ,-10  ,300,  0.1, 1000     ) ==doctest::Approx(5));  
  CHECK(computeNewWTD(6270.33057397168  ,-10  ,300,  0.2, 900      ) ==doctest::Approx(5));  
  CHECK(computeNewWTD(1246.26246881239  ,-3   ,300,  0.5, 500      ) ==doctest::Approx(1));  
  CHECK(computeNewWTD(3892.0177481876   ,-2   ,300,  0.8, 1500     ) ==doctest::Approx(1));  
  CHECK(computeNewWTD(125512.182048359  ,-100 ,300,  0.3, 1000     ) ==doctest::Approx(100));
  CHECK(computeNewWTD(7854.87745892122  ,-10  ,100,  0.3, 1000     ) ==doctest::Approx(5));  
  CHECK(computeNewWTD(7985.04987524957  ,-10  ,1000, 0.3, 1000     ) ==doctest::Approx(5));  
  CHECK(computeNewWTD(6896.36167648567  ,-10  ,10,   0.3, 1000     ) ==doctest::Approx(5));  

  CHECK(computeNewWTD(-487.660600188348  ,-5  ,300, 0.1, 1000      ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-877.789080339026  ,-5  ,300, 0.2, 900       ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-746.262468812392  ,0   ,300, 0.5, 500       ) ==doctest::Approx(-3));
  CHECK(computeNewWTD(-1194.01552781598  ,-1  ,300, 0.8, 1500      ) ==doctest::Approx(-2));
  CHECK(computeNewWTD(-215.318057232321  ,-99 ,300, 0.3, 1000      ) ==doctest::Approx(-100));
  CHECK(computeNewWTD(-1391.76019394264  ,-5  ,100, 0.3, 1000      ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-1488.79363305426  ,-5  ,1000,0.3, 1000      ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-715.953655623573  ,-5  ,10,  0.3, 1000      ) ==doctest::Approx(-10));

  CHECK(computeNewWTD(-4000              ,5   ,300, 0.1, 1000      ) ==doctest::Approx(1));
  CHECK(computeNewWTD(-4500              ,5   ,300, 0.2, 900       ) ==doctest::Approx(0));
  CHECK(computeNewWTD(-5000              ,20  ,300, 0.5, 500       ) ==doctest::Approx(10));
  CHECK(computeNewWTD(-15000             ,10  ,300, 0.8, 1500      ) ==doctest::Approx(0));

  CHECK(computeNewWTD(-5983.51698553982  ,5   ,300,  0.1, 1000     ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-6270.33057397168  ,5   ,300,  0.2, 900      ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-1246.26246881239  ,1   ,300,  0.5, 500      ) ==doctest::Approx(-3));
  CHECK(computeNewWTD(-3892.0177481876   ,1   ,300,  0.8, 1500     ) ==doctest::Approx(-2));
  CHECK(computeNewWTD(-125512.182048359  ,100 ,300,  0.3, 1000     ) ==doctest::Approx(-100));
  CHECK(computeNewWTD(-7854.87745892122  ,5   ,100,  0.3, 1000     ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-7985.04987524957  ,5   ,1000, 0.3, 1000     ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-6896.36167648567  ,5   ,10,   0.3, 1000     ) ==doctest::Approx(-10));

}


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
