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

typedef double* d2d_pointer;
typedef float*  f2d_pointer;
typedef double* d1d_pointer;
typedef float*  f1d_pointer;

struct FanDarcyPack {
  int         width;
  double      cellsize_n_s_metres;
  d1d_pointer cell_area;
  d1d_pointer cellsize_e_w_metres;
  f2d_pointer fdepth;
  f2d_pointer porosity;
  f2d_pointer topo;
  f2d_pointer transmissivity;
  f2d_pointer wtd;
};



double that_one_equation(
  const double fdepth,
  const double wtd,
  const double volume,
  const double capacity,
  const bool gain
){
  if(gain){
    return
      - fdepth
      * std::log( std::exp(wtd / fdepth) - volume / capacity )
      + wtd;
  }
  else{
    return
      - fdepth
      * std::log( std::exp(wtd / fdepth) + volume / capacity )
      + wtd;
  }
}



double that_equation6_thing(
  const double wtd,
  const double fdepth,
  const double ksat
){
  constexpr float shallow = 1.5;
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



double that_basic_diffusion_thing(
  const double Dmax,
  const double cell_width,
  const double cell_height
){
  constexpr double SAFETY = 2;
  return
      (                 std::pow(cell_width, 2) * std::pow(cell_height, 2)) /
      (SAFETY * Dmax * (std::pow(cell_width, 2) + std::pow(cell_height, 2)));
}



double some_other_equation(
  const double capacity,
  const double fdepth,
  const double lower_edge,
  const double upper_edge
){
  return -capacity * (std::exp(lower_edge / fdepth) - std::exp(upper_edge / fdepth));
}

/**
 * @brief Returns the maximum stable time step with a 2x factor of safety
 * @details Uses a 2D diffusion von Neumann stability analysis using the
 *          "worst case" scenario highest transmissivity, combined with
 *          a porosity-based amplification factor.
 */
double computeMaxStableTimeStep(
  const int x,
  const int y,
  const FanDarcyPack &fdp
){
  const int width = fdp.width;

  // Transmissivity is an effective diffusivity
  // Use the highest transmissivity for this time-step calculation

  const std::array<double,4> Tarray = {
    c2d(fdp.transmissivity, x  , y+1),
    c2d(fdp.transmissivity, x  , y-1),
    c2d(fdp.transmissivity, x-1, y  ),
    c2d(fdp.transmissivity, x+1, y  )
  };

  const auto Dmax = *std::max_element(Tarray.begin(), Tarray.end());
  const auto dt_max_diffusion_basic = that_basic_diffusion_thing(Dmax, fdp.cellsize_e_w_metres[y], fdp.cellsize_n_s_metres);

  // Now let's add in the porosity differences
  // Porosity differences will AMPLIFY the WTD changes.
  // Let us choose a time step that is also based on the LOWEST porosity
  // (highest WTD change per water volume transfer)
  const std::array<float,5> PhiArray = {
    c2d(fdp.porosity, x  , y  ),
    c2d(fdp.porosity, x+1, y  ),
    c2d(fdp.porosity, x-1, y  ),
    c2d(fdp.porosity, x  , y+1),
    c2d(fdp.porosity, x  , y-1)
  };

  const auto PhiMin = std::max(0.2f, *std::min_element(PhiArray.begin(), PhiArray.end()));
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



/**
 * @brief Calculates the change in water volume that occurs between
 * two cells, given the water-table depth flux between the two.
 */
double calculateWaterVolume(
  const float wtd_change,
  const float center_wtd,
  const float neighbour_wtd,
  const int x1,
  const int y1,
  const int x2,
  const int y2,
  const FanDarcyPack &fdp
){
  const int width = fdp.width;
  double volume_change = 0;

  // We have the change in the water table between the two cells.
  // We will use this to calculate the volume that moves between the two,
  // based on the average porosity and average water table depths
  // of the two cells.
  // This volume can then be used to calculate the actual
  // change in water table height in each of the two cells.

  const float mean_wtd      = (             center_wtd + neighbour_wtd          ) / 2.;
  const float mean_porosity = (c2d(fdp.porosity,x1,y1) + c2d(fdp.porosity,x2,y2)) / 2.;
  const float mean_fdepth   = (c2d(fdp.fdepth,  x1,y1) + c2d(fdp.fdepth,  x2,y2)) / 2.;
  const float mean_area     = (      fdp.cell_area[y1] + fdp.cell_area[y2]      ) / 2.;

  const float upper_edge = mean_wtd + wtd_change/2.;
  const float lower_edge = mean_wtd - wtd_change/2.;
  //since I am using averages between the two cells, I have to take
  //water both moving up or moving down into account.
  //So, I treat the average wtd as the middle of where water is
  //with respect to changing porosity.

  const auto mean_capacity = mean_area * mean_porosity * mean_fdepth;

  // First, we check to see if the mean wtd was above the surface.

  if(lower_edge > 0){
    // If wtd stays above 0 once the change has occurred, then no need to worry
    // about porosity.

    // The volume change is just the height change multiplied by the cell's
    // area.
    volume_change = wtd_change * mean_area;
  } else if(lower_edge < 0 && upper_edge > 0){
    // The water table drops below the surface during this iteration, so we
    // need to consider porosity for part of the water.

    //TODO: Can this be a named equation?

    // There is a portion of the water that is above the land surface.
    volume_change = upper_edge * mean_area;
    // The portion that is below the land surface.
    // -= because this comes out as a negative number.
    volume_change -= mean_capacity * std::expm1(lower_edge / mean_fdepth);
  } else{
    // the water table is below the surface to start with, therefore it is
    // below the surface the whole time and we need porosity for all the change.
    volume_change = some_other_equation(mean_capacity, mean_fdepth, lower_edge, upper_edge);
  }

  // So now we have the volume change that occurs in each of the two cells
  // as a positive value  whether it was all above ground, all below ground,
  // or a combination.

  return volume_change;
}


/**
 * @brief Calculates water-table depth change in a cell that receives water,
 * given the change in the corresponding cell that gives water.
 */
double computeNewWTDGain(
  const float volume,
  const float my_wtd,
  const float fdepth,
  const float porosity,
  const float cell_area
){

  //We convert the known change in water volume to a change
  //change in water table depth for this cell.

  //We have a volume, which is positive; we also know that this cell is gaining water.

  const auto gw_storage_cap = cell_area * porosity * fdepth; //Groundwater Storage Capacity

  double change_in_wtd;

  if(my_wtd < 0){
    //either all of the change is below the surface,
    //or it's a combination of above and below.
    //we start off assuming that it will all be below the surface:
    change_in_wtd = that_one_equation(fdepth, my_wtd, volume, gw_storage_cap, true);

    //TODO(kerry): Find names for this conditions
    const bool cond1 = my_wtd + change_in_wtd > 0;
    const bool cond2 = std::exp(my_wtd / fdepth) < volume / gw_storage_cap; //TODO: Should this be expm1, like below?
    if(cond1 || cond2){
      //TODO: Should this be that_on_equation?

      //it is gaining enough water that some will be above the surface.
      //we want to calculate how much of the water is used up
      //in the ground, i.e. the portion between the starting wtd and 0.
      const double GW_portion = -gw_storage_cap * std::expm1(my_wtd/fdepth); //TODO: Should this be const? It's subtracted below
      //this is the volume of water used up in filling in the ground.
      //the volume minus GW_portion is left over to fill surface water.
      change_in_wtd = (volume - GW_portion) / cell_area - my_wtd;
      //-my_wtd because this was a negative number and we are getting the total change in wtd.
    }
  } else {
    //all of the change is above the surface.
    change_in_wtd = volume / cell_area;
  }

  return change_in_wtd;
}


/**
 * @brief Calculates water-table depth change in a cell that gives water,
 * given the change in the corresponding cell that receives water.
 */
double computeNewWTDLoss(
  const float volume,
  const float my_wtd,
  const float fdepth,
  const float porosity,
  const float cell_area
){
  //We convert the known change in water volume to a change
  //change in water table depth for this cell.

  //We have a volume, which is positive; we also know that this cell is losing water.

  const auto gw_storage_cap = cell_area * porosity * fdepth; //Groundwater Storage Capacity

  double change_in_wtd = -volume / cell_area;

  if((my_wtd > 0) && (my_wtd + change_in_wtd < 0)){
    //at least some of the change is above the surface.
    //how much of the volume is used up above the surface?
    const double SW_portion = my_wtd * cell_area;
    //this is the volume used in water above the surface.
    //The total volume minus this is left over to decrease groundwater.
    change_in_wtd = that_one_equation(fdepth, my_wtd, volume-SW_portion, gw_storage_cap, false);
    //the cell is losing water, and it's losing enough
    //that some of the change is below the surface and
    //so we need to take porosity into account.
  } else if(my_wtd < 0){
    //Since it's losing water, all of the change is below the surface.
    change_in_wtd = that_one_equation(fdepth, my_wtd, volume, gw_storage_cap, false);
  }

  return change_in_wtd;
}



struct double2 {
  double x = 0;
  double y = 0;
  double2() = default;
  double2(double x, double y) : x(x), y(y) {}
  double2& operator+=(const double2 &o) {
    x+=o.x;
    y+=o.y;
    return *this;
  }
};



double2 GainLoss(
  const int x,  //Focal cell
  const int y,
  const int nx, //Neighbour cell
  const int ny,
  const double wtd_change,
  const double volume,
  const double local_wtd_center,
  const double local_wtd_neighbour,
  const FanDarcyPack &fdp
){
  const int width = fdp.width;

  double local_wtd_n_change;
  double mycell_change;
  if(wtd_change > 0){  // The current cell is giving water to the neighbour.
    local_wtd_n_change = computeNewWTDLoss(volume, local_wtd_neighbour,
                            c2d(fdp.fdepth,    nx, ny),
                            c2d(fdp.porosity,  nx, ny),
                            c2d(fdp.cell_area, nx, ny)
                          );
    mycell_change      = computeNewWTDGain(volume, local_wtd_center,
                            c2d(fdp.fdepth,     x, y ),
                            c2d(fdp.porosity,   x, y ),
                            c2d(fdp.cell_area,  x, y )
                          );
  } else {             // The neighbour is the receiving cell.
    local_wtd_n_change = computeNewWTDGain(volume, local_wtd_neighbour,
                            c2d(fdp.fdepth,    nx, ny),
                            c2d(fdp.porosity,  nx, ny),
                            c2d(fdp.cell_area, nx, ny)
                          );
    mycell_change      = computeNewWTDLoss(volume, local_wtd_center,
                            c2d(fdp.fdepth,    x,  y ),
                            c2d(fdp.porosity,  x,  y ),
                            c2d(fdp.cell_area, x,  y )
                          );
  }

  return double2(mycell_change, local_wtd_n_change);
}



//TODO: Maybe use this
double2 computeWTDchangeWithNeighbour(
  const int x,  //Focal cell coordinates
  const int y,
  const int nx, //Neighbour cell coordinates
  const int ny,
  const double wtd,
  const double wtd_n,
  const double dt,
  const std::array<double,5> &local_wtd,
  const bool north_south,
  const FanDarcyPack &fdp
){
  const auto width = fdp.width;
  double Q;

  const auto headCenter = c2d(fdp.topo, x,  y ) + wtd;
  const auto head_n     = c2d(fdp.topo, nx, ny) + wtd_n;

  const auto transmissivity = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, nx, ny)) / 2.;

  if(north_south)
    Q = transmissivity * (head_n - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  else
    Q = transmissivity * (head_n - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;


  const auto wtd_change = Q * dt / fdp.cell_area[ny];

  const auto volume = calculateWaterVolume(std::abs(wtd_change), wtd, wtd_n, x, y, nx, ny, fdp);

  if(volume > 0){
    return GainLoss(x, y, nx  , ny, wtd_change, volume, wtd, wtd_n, fdp);
  } else {
    return {};
  }
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
  const double transmissivityE = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x+1, y  )) / 2.;
  const double transmissivityW = (c2d(fdp.transmissivity, x, y) + c2d(fdp.transmissivity, x-1, y  )) / 2.;

  // Then, compute the discharges
  const double QN = transmissivityN * (headN - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QS = transmissivityS * (headS - headCenter) / fdp.cellsize_n_s_metres    * fdp.cellsize_e_w_metres[y];
  const double QE = transmissivityE * (headE - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;
  const double QW = transmissivityW * (headW - headCenter) / fdp.cellsize_e_w_metres[y] * fdp.cellsize_n_s_metres;

  // Update water-table depth, but only in the "internal" variables to
  // handle internal time stepping to maintain stability.
  // dH = sum(discharges) times time step, divided by cell area,
  //      divided by porosity.
  const double wtd_change_N = QN * dt / fdp.cell_area[y+1];
  const double wtd_change_S = QS * dt / fdp.cell_area[y-1];
  const double wtd_change_W = QW * dt / fdp.cell_area[y  ];
  const double wtd_change_E = QE * dt / fdp.cell_area[y  ];

  //We use the calculated wtd change fluxes to calculate a volume
  //change between the 2 cells. This is based on average
  //porosity, e-folding depth, wtd, and cell area of the 2 cells.

  const double volume_N = calculateWaterVolume(std::abs(wtd_change_N), local_wtd[0], local_wtd[1], x, y, x, y+1, fdp);
  const double volume_S = calculateWaterVolume(std::abs(wtd_change_S), local_wtd[0], local_wtd[2], x, y, x, y-1, fdp);
  const double volume_E = calculateWaterVolume(std::abs(wtd_change_E), local_wtd[0], local_wtd[3], x, y, x+1, y, fdp);
  const double volume_W = calculateWaterVolume(std::abs(wtd_change_W), local_wtd[0], local_wtd[4], x, y, x-1, y, fdp);

  //we now use these volumes to compute the actual changes in
  //water table depths in the target cell and each of the neighbouring cells:

  //The current cell is receiving water from the North, so (x,y+1) is the giving cell.
  //Target cell is the receiving cell.
  double my_cell_change = 0;
  if(volume_N > 0){
    const auto change = GainLoss(x, y, x  , y+1, wtd_change_N, volume_N, local_wtd[0], local_wtd[1], fdp);
    my_cell_change += change.x;
    local_wtd[1]   += change.y;
  }
  if(volume_S > 0){
    const auto change = GainLoss(x, y, x  , y-1, wtd_change_S, volume_S, local_wtd[0], local_wtd[2], fdp);
    my_cell_change += change.x;
    local_wtd[2]   += change.y;
  }
  if(volume_E > 0){
    const auto change = GainLoss(x, y, x+1, y  , wtd_change_E, volume_E, local_wtd[0], local_wtd[3], fdp);
    my_cell_change += change.x;
    local_wtd[3]   += change.y;
  }
  if(volume_W > 0){
    const auto change = GainLoss(x, y, x-1, y  , wtd_change_W, volume_W, local_wtd[0], local_wtd[4], fdp);
    my_cell_change += change.x;
    local_wtd[4]   += change.y;
  }

  // Now we have the height changes that will take place in the target cell
  // and each of the four neighbours.
  local_wtd[0] += my_cell_change;
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

void update(const Parameters &params, ArrayPack &arp){
  FanDarcyPack fdp;
  fdp.cell_area           = arp.cell_area.data();
  fdp.cellsize_e_w_metres = arp.cellsize_e_w_metres.data();
  fdp.cellsize_n_s_metres = params.cellsize_n_s_metres;
  fdp.fdepth              = arp.fdepth.data();
  fdp.porosity            = arp.porosity.data();
  fdp.topo                = arp.topo.data();
  fdp.transmissivity      = arp.transmissivity.data();
  fdp.width               = arp.fdepth.width();
  fdp.wtd                 = arp.wtd.data();

  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    arp.transmissivity(x,y) = that_equation6_thing(
      arp.wtd(x,y),
      arp.fdepth(x,y),
      arp.ksat(x,y)
    );
  }

  // Updates water-table depth grid over one time step
  #pragma omp parallel for collapse(2) default(none) shared(arp,params,fdp)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    // Skip ocean cells
    if(arp.land_mask(x,y) == 1){
      // Otherwise, update the water-table depth change array
      arp.wtd_changed(x,y) = updateCell(x, y, params.deltat, fdp);
    }
  }


  // Once all the changes are known, update the WTD everywhere with the
  // difference array
  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
  // Skip ocean cells
    if(arp.land_mask(x,y) == 1){
      // Update the whole wtd array at once.
      // This is the new water table after groundwater has moved
      // for delta_t seconds.
      arp.wtd(x,y) = arp.wtd_changed(x,y);;
    }
  }

}

}



TEST_CASE("that_one_equation"){
  CHECK(FanDarcyGroundwater::that_one_equation(1,1,1,2,true)==doctest::Approx(0.653983278399143));
  CHECK(FanDarcyGroundwater::that_one_equation(1,4,3,4,true)==8);
  CHECK(FanDarcyGroundwater::that_one_equation(1,2,3,6,true)==9);
}

TEST_CASE("that_equation6_thing"){
  CHECK(FanDarcyGroundwater::that_equation6_thing(1,3,6)==9);
}
