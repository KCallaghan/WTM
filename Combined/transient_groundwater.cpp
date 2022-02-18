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
   // const double Hmax,

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

  // Use the smallest cellsize of the cells examined for the time-step calculation
  const std::array<double,3> sizeArray = {
    fdp.cellsize_e_w_metres[y],
    fdp.cellsize_e_w_metres[y+1],
    fdp.cellsize_e_w_metres[y-1]
  };

  const auto Dmax = *std::max_element(Tarray.begin(), Tarray.end());
  const auto sizeMin = *std::min_element(sizeArray.begin(), sizeArray.end());

  const auto dt_max_diffusion_basic = vonNeumannStability(Dmax, sizeMin, fdp.cellsize_n_s_metres);


  // Now let's add in the porosity differences
  // Porosity differences will AMPLIFY the WTD changes.
  // Let us choose a time step that is also based on the LOWEST porosity
  // (highest WTD change per water volume transfer)
  auto PhiMin =             c2d(fdp.porosity, x  , y  );
  PhiMin = std::min(PhiMin, c2d(fdp.porosity, x  , y+1  ));
  PhiMin = std::min(PhiMin, c2d(fdp.porosity, x  , y-1  ));
  PhiMin = std::min(PhiMin, c2d(fdp.porosity, x+1  , y  ));
  PhiMin = std::min(PhiMin, c2d(fdp.porosity, x-1  , y  ));
  PhiMin = std::max(PhiMin, 0.2f); //Yes, we do want max - porosity should be a minumum of 0.2 in all cells.

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
  const double porosity,
  const double cell_area
){
  //Since we are using a vertically constant porosity, we can just convert the
  //volume to a height of water right away.
  //If above ground, this will be the actual change in height;
  //If below ground, dividing this by porosity will give the height change.
  const double porosity_one_height_change = volume_change/cell_area;
  double final_wtd;

  //either the cell starts with wtd above the surface, or with it below the surface.
  //If above the surface:
  if(initial_wtd >= 0){
    //Either the result after moving water still has water above the surface,
    //because it was gaining water, or because it was losing water but it
    //lost less than the total surface water, or
    //the results after moving water has the water table below the surface.

    //Now, let's compare the initial wtd + the height change to the land surface:
    //If positive, it is still above the surface,
    //if negative, it is below the surface.
    final_wtd = initial_wtd + porosity_one_height_change;  //if positive, this is the final answer.
    if(final_wtd < 0){
      //the water table will now be below the surface, so we have to take porosity into account.
      //final_wtd currently represents where the water table would be if porosity were 1.
      //all we have to do is divide it by porosity to get the actual water table.
      final_wtd /= porosity;
    }

  }
  else{
    //Otherwise, if the cell starts with wtd below the surface:
    //what height change would be enough to switch to above the surface?
    //let's represent the dry height below-ground as it would have been with porosity = 1:
    double porosity_one_below_ground_height = -initial_wtd * porosity;
    if(porosity_one_height_change > porosity_one_below_ground_height){
      //the final wtd will be greater than 0. Subtract the amount used up below ground
      //to get the amount left above ground.
      final_wtd = porosity_one_height_change - porosity_one_below_ground_height;
    }
    else{
      //the final wtd is still below ground; just divide the height change by
      //porosity to get the below-ground height change.
      final_wtd = initial_wtd + (porosity_one_height_change/porosity);
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

  const double transmissivityN = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x  , y+1)) / 2.;
  const double transmissivityS = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x  , y-1)) / 2.;
  const double transmissivityE = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x+1  , y)) / 2.;
  const double transmissivityW = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x-1  , y)) / 2.;

  // Then, compute the discharges
  // Define Q such that + if center cell gaining, - if center cell losing
  const double QN = transmissivityN * (headN - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QS = transmissivityS * (headS - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QE = transmissivityE * (headE - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;
  const double QW = transmissivityW * (headW - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;

  // Sum discharges and divide by the time step to compute the water volume
  // added to the center cell
  const double dVolume = (QN + QS + QW + QE) * dt;
if(!(local_wtd[0]<0 || local_wtd[0]>=0)){
  std::cout<<"cell was "<<x<<" "<<y<<std::endl;
  std::cout<<"wtd in the cell was "<<local_wtd[0]<<std::endl;
  std::cout<<"transmissivityN "<<transmissivityN<<" transmissivityS "<<transmissivityS<<" transmissivityE "<<transmissivityE<<" transmissivityW "<<transmissivityW<<std::endl;
  std::cout<<"headCenter "<<headCenter<<" headN "<<headN<<" headS "<<headS<<" headE "<<headE<<" headW "<<headW<<std::endl;
  std::cout<<"QN "<<QN<<" QS "<<QS<<" QE "<<QE<<" QW "<<QW<<" dVolume "<<dVolume<<std::endl;
}
  // Update the cell's WTD
  local_wtd[0] = computeNewWTD( dVolume, local_wtd[0], c2d(fdp.porosity, x, y), fdp.cell_area[y] );

if(!(local_wtd[0]<0 || local_wtd[0]>=0)){
  std::cout<<"cell was "<<x<<" "<<y<<std::endl;
  std::cout<<"topo was "<<c2d(fdp.topo, x  , y  )<<std::endl;
  std::cout<<"wtd in the cell now is "<<local_wtd[0]<<std::endl;
  std::cout<<"transmissivityN "<<transmissivityN<<" transmissivityS "<<transmissivityS<<" transmissivityE "<<transmissivityE<<" transmissivityW "<<transmissivityW<<std::endl;
  std::cout<<"headCenter "<<headCenter<<" headN "<<headN<<" headS "<<headS<<" headE "<<headE<<" headW "<<headW<<std::endl;
  std::cout<<"QN "<<QN<<" QS "<<QS<<" QE "<<QE<<" QW "<<QW<<" dVolume "<<dVolume<<std::endl;

}
  // For local dynamic time stepping (consider switching to global later),
  // update the neighboring cell WTDs
  local_wtd[1] = computeNewWTD( -QN*dt, local_wtd[1], c2d(fdp.porosity, x, y+1), fdp.cell_area[y+1] );
  local_wtd[2] = computeNewWTD( -QS*dt, local_wtd[2], c2d(fdp.porosity, x, y-1), fdp.cell_area[y-1] );
  local_wtd[3] = computeNewWTD( -QE*dt, local_wtd[3], c2d(fdp.porosity, x+1, y), fdp.cell_area[y] );
  local_wtd[4] = computeNewWTD( -QW*dt, local_wtd[4], c2d(fdp.porosity, x-1, y), fdp.cell_area[y] );


}



/**
 * @brief Updates the wtd_depth_total array at a cell(x,y) using the
 * pre-set time step and dynamic time stepping within this as needed.
 */
double updateCell(
  const int x,
  const int y,
  double time_remaining,
  const FanDarcyPack &fdp,
  Parameters &params
){
  const auto width = fdp.width;

  // Skip ocean cells
  if(c2d(fdp.land_mask,x,y) != 1 || c2d(fdp.ice_mask,x,y) != 0)
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

  if( c2d(fdp.land_mask,x+1,y) == 0)
    params.total_loss_to_ocean += local_wtd[3]*fdp.cell_area[y];
  if(c2d(fdp.land_mask,x-1,y) == 0)
    params.total_loss_to_ocean += local_wtd[4]*fdp.cell_area[y];
  if(c2d(fdp.land_mask,x,y+1) == 0)
    params.total_loss_to_ocean += local_wtd[1]*fdp.cell_area[y+1];
  if(c2d(fdp.land_mask,x,y-1) == 0)
    params.total_loss_to_ocean += local_wtd[2]*fdp.cell_area[y-1];

  // When exiting loop, the local_wtd[0] variable holds the final water-table depth
  return local_wtd[0];
}



//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void UpdateCPU(Parameters &params, ArrayPack &arp){
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


  // Updates water-table depth grid over one time step
  #pragma omp parallel for collapse(2) default(none) shared(arp,params,fdp)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    arp.wtd_changed(x,y) = updateCell(x, y, params.deltat, fdp,params);
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



void update(Parameters &params, ArrayPack &arp){
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
  CHECK(computeNewWTD(6000,   5, 0.4, 1000 ) ==doctest::Approx(11));
  CHECK(computeNewWTD(6000,   1, 0.4, 1000 ) ==doctest::Approx(7));
  CHECK(computeNewWTD(6000,   1, 0.4, 10000) ==doctest::Approx(1.6));
  CHECK(computeNewWTD(6000,   1, 0.4, 100  ) ==doctest::Approx(61));
  CHECK(computeNewWTD(6000,   1, 0.2, 1000 ) ==doctest::Approx(7));
  CHECK(computeNewWTD(6000,   1, 0.5, 1000 ) ==doctest::Approx(7));
  CHECK(computeNewWTD(6000,   1, 0.8, 1000 ) ==doctest::Approx(7));
  CHECK(computeNewWTD(100000, 1, 0.4, 1000 ) ==doctest::Approx(101));
  CHECK(computeNewWTD(1000,   1, 0.4, 1000 ) ==doctest::Approx(2));
  CHECK(computeNewWTD(  10,   1, 0.4, 1000 ) ==doctest::Approx(1.01));

  CHECK(computeNewWTD(-6000,   5, 0.4, 1000 ) ==doctest::Approx(-2.5));
  CHECK(computeNewWTD(-6000,   1, 0.4, 1000 ) ==doctest::Approx(-12.5));
  CHECK(computeNewWTD(-6000,   1, 0.4, 10000) ==doctest::Approx(0.4));
  CHECK(computeNewWTD(-6000,   1, 0.4, 100  ) ==doctest::Approx(-147.5));
  CHECK(computeNewWTD(-6000,   1, 0.2, 1000 ) ==doctest::Approx(-25));
  CHECK(computeNewWTD(-6000,   1, 0.5, 1000 ) ==doctest::Approx(-10));
  CHECK(computeNewWTD(-6000,   1, 0.8, 1000 ) ==doctest::Approx(-6.25));
  CHECK(computeNewWTD(-100000, 1, 0.4, 1000 ) ==doctest::Approx(-247.5));
  CHECK(computeNewWTD(-1000,   1, 0.4, 1000 ) ==doctest::Approx(0));
  CHECK(computeNewWTD(  -10,   1, 0.4, 1000 ) ==doctest::Approx(0.99));

  CHECK(computeNewWTD(6000,   -5, 0.4, 1000 ) ==doctest::Approx(4));
  CHECK(computeNewWTD(6000,   -1, 0.4, 1000 ) ==doctest::Approx(5.6));
  CHECK(computeNewWTD(6000,   -1, 0.4, 10000) ==doctest::Approx(0.2));
  CHECK(computeNewWTD(6000,   -1, 0.4, 100  ) ==doctest::Approx(59.6));
  CHECK(computeNewWTD(6000,   -1, 0.2, 1000 ) ==doctest::Approx(5.8));
  CHECK(computeNewWTD(6000,   -1, 0.5, 1000 ) ==doctest::Approx(5.5));
  CHECK(computeNewWTD(6000,   -1, 0.8, 1000 ) ==doctest::Approx(5.2));
  CHECK(computeNewWTD(100000, -1, 0.4, 1000 ) ==doctest::Approx(99.6));
  CHECK(computeNewWTD(1000,   -1, 0.4, 1000 ) ==doctest::Approx(0.6));
  CHECK(computeNewWTD(  10,   -1, 0.4, 1000 ) ==doctest::Approx(-0.975));

  CHECK(computeNewWTD(-6000,   -5, 0.4, 1000 ) ==doctest::Approx(-20));
  CHECK(computeNewWTD(-6000,   -1, 0.4, 1000 ) ==doctest::Approx(-16));
  CHECK(computeNewWTD(-6000,   -1, 0.4, 10000) ==doctest::Approx(-2.5));
  CHECK(computeNewWTD(-6000,   -1, 0.4, 100  ) ==doctest::Approx(-151));
  CHECK(computeNewWTD(-6000,   -1, 0.2, 1000 ) ==doctest::Approx(-31));
  CHECK(computeNewWTD(-6000,   -1, 0.5, 1000 ) ==doctest::Approx(-13));
  CHECK(computeNewWTD(-6000,   -1, 0.8, 1000 ) ==doctest::Approx(-8.5));
  CHECK(computeNewWTD(-100000, -1, 0.4, 1000 ) ==doctest::Approx(-251));
  CHECK(computeNewWTD(-1000,   -1, 0.4, 1000 ) ==doctest::Approx(-3.5));
  CHECK(computeNewWTD(  -10,   -1, 0.4, 1000 ) ==doctest::Approx(-1.025));

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
