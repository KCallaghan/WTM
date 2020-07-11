#include <algorithm>
#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"

const double FP_ERROR = 1e-4;

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {


/**
 * @brief Returns the maximum stable time step with a 2x factor of safety
 * @details Uses a 2D diffusion von Neumann stability analysis using the
 *          "worst case" scenario highest transmissivity, combined with
 *          a porosity-based amplification factor.
 */
double computeMaxStableTimeStep(const Parameters &params,
                                                     ArrayPack &arp,
                                                     uint32_t x,
                                                     uint32_t y){
  // Transmissivity is an effective diffusivity
  // Use the highest transmissivity for this time-step calculation
  double dt_max_diffusion_basic;
  double dt_max_diffusion_withPorosity;

  std::array<double,4> Tarray = { arp.transmissivity(x,y+1), arp.transmissivity(x,y-1),
                           arp.transmissivity(x-1,y), arp.transmissivity(x+1,y) };
  const double Dmax = *std::max_element(Tarray.begin(), Tarray.end());
  dt_max_diffusion_basic = ( pow(arp.cellsize_e_w_metres[y], 2)
                                * pow(params.cellsize_n_s_metres, 2) ) \
                                / ( 2 * Dmax \
                                  * ( pow(arp.cellsize_e_w_metres[y], 2) \
                                    + pow(params.cellsize_n_s_metres, 2) ) );

  // Now let's add in the porosity differences
  // Porosity differences will AMPLIFY the WTD changes.
  // Let us choose a time step that is also based on the LOWEST porosity
  // (highest WTD change per water volume transfer)
  std::array<float,5> PhiArray = { arp.porosity(x, y),
                                   arp.porosity(x+1, y), arp.porosity(x-1, y),
                                   arp.porosity(x, y+1), arp.porosity(x, y-1) };
  float PhiMin = *std::min_element(PhiArray.begin(), PhiArray.end());
  if(PhiMin < 0.2)
      PhiMin = 0.2;
  // Porosity is a linear amplifier of WTD change, and it amplifies change
  // in both the giving and receiving cells.
  // Amplification goes as 1 / phi.
  // We need to consider this going up and down, so the amplification becomes
  // 2/Phi, and therefore we multiply our stable dt by Phi/2
  dt_max_diffusion_withPorosity = dt_max_diffusion_basic * PhiMin / 2.;
  // In order to avoid operating at the very maximum time step possible,
  // we apply a factor of safety of 2

  return dt_max_diffusion_withPorosity/2.;
}


/**
 * @brief Calculates the change in water volume that occurs between
 * two cells, given the water-table depth flux between the two.
 */
double calculateWaterVolume(const float wtd_change,
                                                 const float center_wtd,
                                                 const float neighbour_wtd,
                                                 const int x1,
                                                 const int y1,
                                                 const int x2,
                                                 const int y2,
                                                 const ArrayPack &arp){


  double volume_change = 0.;

   // We have the change in the water table between the two cells.
   // We will use this to calculate the volume that moves between the two,
   // based on the average porosity and average water table depths
   // of the two cells.
   // This volume can then be used to calculate the actual
   // change in water table height in each of the two cells.

  float mean_wtd = (center_wtd + neighbour_wtd) / 2.;
  float mean_porosity = (arp.porosity(x1,y1) + arp.porosity(x2,y2)) / 2.;
  float mean_fdepth = (arp.fdepth(x1,y1) + arp.fdepth(x2,y2)) / 2.;
  float mean_area = (arp.cell_area[y1] + arp.cell_area[y2]) / 2.;

  float upper_edge = mean_wtd + wtd_change/2.;
  float lower_edge = mean_wtd - wtd_change/2.;
  //since I am using averages between the two cells, I have to take
  //water both moving up or moving down into account.
  //So, I treat the average wtd as the middle of where water is
  //with respect to changing porosity.


  // First, we check to see if the mean wtd was above the surface.

  // If wtd stays above 0 once the change has occurred, then no need to worry
  // about porosity.
  if(lower_edge > 0){
    // The volume change is just the height change multiplied by the cell's
    // area.
    volume_change = wtd_change * mean_area;
    // The water table drops below the surface during this iteration, so we
    // need to consider porosity for part of the water.
  }
  else if(lower_edge < 0 && upper_edge > 0){
  // There is a portion of the water that is above the land surface.
    volume_change = upper_edge * mean_area;
    // The portion that is below the land surface.
    // -= because this comes out as a negative number.
    volume_change -= mean_area * mean_porosity * mean_fdepth * \
                        ( exp( (lower_edge) / mean_fdepth ) \
                          - 1);

  }
  // the water table is below the surface to start with, therefore it is
  // below the surface the whole time and we need porosity for all the change.
  else{
    volume_change = -mean_area * mean_porosity * mean_fdepth * \
                      ( exp(lower_edge / mean_fdepth) - \
                        exp(upper_edge / mean_fdepth) );
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
double computeNewWTDGain(const float volume,
                                              const float my_wtd,
                                              const int x,
                                              const int y,
                                              const ArrayPack &arp){

  //We convert the known change in water volume to a change
  //change in water table depth for this cell.

  //We have a volume, which is positive; we also know that this cell is gaining water.

  //at least some of the change is above the surface.
  double change_in_wtd = volume / arp.cell_area[y];

  if(my_wtd < 0){
    //either all of the change is below the surface,
    //or it's a combination of above and below.
    //we start off assuming that it will all be below the surface:
    change_in_wtd = -arp.fdepth(x,y) * log( exp(my_wtd / arp.fdepth(x,y)) \
                    - volume / (arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)) ) + my_wtd;

    if((my_wtd + change_in_wtd > 0) || exp(my_wtd / arp.fdepth(x,y)) < \
        (volume/(arp.cell_area[y] * arp.porosity(x,y)*arp.fdepth(x,y)))  ){
      //it is gaining enough water that some will be above the surface.
      //we want to calculate how much of the water is used up
      //in the ground, i.e. the portion between the starting wtd and 0.
      double GW_portion = -arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y) \
                                  * ( exp( my_wtd /arp.fdepth(x,y) ) - 1);
      //this is the volume of water used up in filling in the ground.
      //the volume minus GW_portion is left over to fill surface water.
      change_in_wtd = ( (volume - GW_portion) / arp.cell_area[y] ) - my_wtd;
      //-my_wtd because this was a negative number and we are getting the total change in wtd.
    }
  }

  return change_in_wtd;
}


/**
 * @brief Calculates water-table depth change in a cell that gives water,
 * given the change in the corresponding cell that receives water.
 */
double computeNewWTDLoss(const float volume,
                                              const float my_wtd,
                                              const int x,
                                              const int y,
                                              const ArrayPack &arp){

    //We convert the known change in water volume to a change
    //change in water table depth for this cell.

    //We have a volume, which is positive; we also know that this cell is losing water.

  double change_in_wtd = -volume / arp.cell_area[y];

  if((my_wtd > 0) && (my_wtd + change_in_wtd < 0)){
    //at least some of the change is above the surface.
    //how much of the volume is used up above the surface?
    double SW_portion = my_wtd * arp.cell_area[y];
    //this is the volume used in water above the surface.
    //The total volume minus this is left over to decrease groundwater.
    change_in_wtd = -arp.fdepth(x,y) * log( exp(my_wtd / arp.fdepth(x,y)) \
                    + (volume - SW_portion) / (arp.cell_area[y] * arp.porosity(x,y) \
                    * arp.fdepth(x,y)) ) + my_wtd;

    //the cell is losing water, and it's losing enough
    //that some of the change is below the surface and
    //so we need to take porosity into account.

    }
    else if(my_wtd < 0){
      //Since it's losing water, all of the change is below the surface.
      change_in_wtd = -arp.fdepth(x,y) * log( exp(my_wtd / arp.fdepth(x,y)) \
                      + volume / (arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y))) + my_wtd;
    }

  return change_in_wtd;
}


/**
 * @brief Calculates water-table depth change at a cell (and associated
 *        surrounding cells) and updates the class variables associated
 *        with these.
 */
void computeWTDchangeAtCell( const Parameters &params,
                                                  ArrayPack &arp,
                                                  int32_t x, int32_t y,
                                                  double dt, std::array<double,5> &local_wtd){
  // Update WTD change

  // We do this instead of using a staggered grid to approx. double CPU time
  // in exchange for using less memory.

  // First, compute elevation head at center cell and all neighbours
  // This equals topography plus the water table depth
  // (positive upwards; negative if water table is below Earth surface)
  double headCenter = arp.topo(x,y)   + local_wtd[0];
  double headN      = arp.topo(x,y+1) + local_wtd[1];
  double headS      = arp.topo(x,y-1) + local_wtd[2];
  double headE      = arp.topo(x+1,y) + local_wtd[3];
  double headW      = arp.topo(x-1,y) + local_wtd[4];

  double mycell_change = 0.;

  // Then, compute the discharges
  double QN = ((arp.transmissivity(x,y) + arp.transmissivity(x,y+1)) / 2.) \
                * (headN - headCenter) / params.cellsize_n_s_metres * arp.cellsize_e_w_metres[y];
  double QS = ((arp.transmissivity(x,y) + arp.transmissivity(x,y-1)) / 2.) \
                * (headS - headCenter) / params.cellsize_n_s_metres * arp.cellsize_e_w_metres[y];
  double QE = ((arp.transmissivity(x,y) + arp.transmissivity(x+1,y)) / 2.) \
                * (headE - headCenter) / arp.cellsize_e_w_metres[y] * params.cellsize_n_s_metres;
  double QW = ((arp.transmissivity(x,y) + arp.transmissivity(x-1,y)) / 2.) \
                * (headW - headCenter) / arp.cellsize_e_w_metres[y] * params.cellsize_n_s_metres;

  // Update water-table depth, but only in the "internal" variables to
  // handle internal time stepping to maintain stability.
  // dH = sum(discharges) times time step, divided by cell area,
  //      divided by porosity.
  double wtd_change_N = QN * dt / ( arp.cell_area[y+1]);
  double wtd_change_S = QS * dt / ( arp.cell_area[y-1]);
  double wtd_change_W = QW * dt / ( arp.cell_area[y] );
  double wtd_change_E = QE * dt / ( arp.cell_area[y] );

  //We use the calculated wtd change fluxes to calculate a volume
  //change between the 2 cells. This is based on average
  //porosity, e-folding depth, wtd, and cell area of the 2 cells.

  double volume_N = calculateWaterVolume(fabs(wtd_change_N),
     local_wtd[0], local_wtd[1], x, y, x, y+1, arp);

  double volume_S = calculateWaterVolume(fabs(wtd_change_S),
     local_wtd[0], local_wtd[2], x, y, x, y-1, arp);

  double volume_E = calculateWaterVolume(fabs(wtd_change_E),
     local_wtd[0], local_wtd[3], x, y, x+1, y, arp);

  double volume_W = calculateWaterVolume(fabs(wtd_change_W),
   	 local_wtd[0], local_wtd[4], x, y, x-1, y, arp);


  //we now use these volumes to compute the actual changes in
  //water table depths in the target cell and each of the neighbouring cells:

  // The current cell is receiving water from the North, so (x,y+1) is the giving cell.
  // Target cell is the receiving cell.
  if(volume_N > 0){
    if(wtd_change_N > 0){
      local_wtd[1]   += computeNewWTDLoss( volume_N, local_wtd[1], x, y+1,arp);
      mycell_change  += computeNewWTDGain( volume_N, local_wtd[0], x, y,  arp);
    }
    // The current cell is giving water to the North.
    // The North is the receiving cell.
    else{
      local_wtd[1]   += computeNewWTDGain( volume_N, local_wtd[1], x, y+1, arp);
      mycell_change  += computeNewWTDLoss( volume_N, local_wtd[0], x, y,   arp);
    }
  }

  if(volume_S > 0){
    if(wtd_change_S > 0){
      local_wtd[2]  += computeNewWTDLoss( volume_S, local_wtd[2], x, y-1, arp);
      mycell_change += computeNewWTDGain( volume_S, local_wtd[0], x, y,   arp);
    }
    else {
      local_wtd[2]   += computeNewWTDGain( volume_S, local_wtd[2], x, y-1, arp);
      mycell_change  += computeNewWTDLoss( volume_S, local_wtd[0], x, y,   arp);
    }
  }

  if(volume_E > 0){
    if(wtd_change_E > 0){
      local_wtd[3]   += computeNewWTDLoss( volume_E, local_wtd[3], x+1, y, arp);
      mycell_change  += computeNewWTDGain( volume_E, local_wtd[0], x,   y, arp);
    }
    else {
      local_wtd[3]   += computeNewWTDGain( volume_E, local_wtd[3], x+1, y, arp);
      mycell_change  += computeNewWTDLoss( volume_E, local_wtd[0], x,   y, arp);
    }
  }

  if(volume_W > 0){
    if(wtd_change_W > 0){
      local_wtd[4]   += computeNewWTDLoss( volume_W, local_wtd[4], x-1, y, arp);
      mycell_change  += computeNewWTDGain( volume_W, local_wtd[0], x,   y, arp);
    }
    else{
      local_wtd[4]   += computeNewWTDGain( volume_W, local_wtd[4], x-1, y, arp);
      mycell_change  += computeNewWTDLoss( volume_W, local_wtd[0], x,   y, arp);
    }
  }
  // Now we have the height changes that will take place in the target cell
  // and each of the four neighbours.
  local_wtd[0] += mycell_change;
}


/**
 * @brief Updates the wtd_depth_total array at a cell(x,y) using the
 * pre-set time step and dynamic time stepping within this as needed.
 */
void updateCell( const Parameters &params, ArrayPack &arp,
                                      uint32_t x, uint32_t y ){

  using namespace std::this_thread;     // sleep_for, sleep_until
  using namespace std::chrono_literals; // ns, us, ms, s, h, etc.

  // Runs functions to compute time steps and update WTD for the center cell
  // and its neighbours until the outer time step has been completed

  // Initialize variables for dynamic time stepping
  double time_remaining = params.deltat;
  double dt_inner;

  // Initial water-table depths, prior to updating
  double wtdCenter_initial = arp.wtd(x,y);

  std::array<double, 5> local_wtd = {arp.wtd(x,y),arp.wtd(x,y+1),arp.wtd(x,y-1),arp.wtd(x+1,y),arp.wtd(x-1,y)};

  // Update water-table depths using dynamic time stepping
  while (time_remaining > 0){

    double max_stable_time_step = computeMaxStableTimeStep( params, arp, x, y);
    // Choose the inner-loop time step
    if(time_remaining <= max_stable_time_step){
      dt_inner = time_remaining;
    }
    else{
      dt_inner = max_stable_time_step;
    }

    computeWTDchangeAtCell(params, arp, x, y, dt_inner,local_wtd);
    time_remaining -= dt_inner;
  }

    // When exiting loop, the local_wtd[0] variable holds the final water-table depth
    arp.wtd_changed(x,y) = local_wtd[0];
}



//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void update(const Parameters &params, ArrayPack &arp){

  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){
    if(arp.fdepth(x,y)>0){
      // Equation S6 from the Fan paper
      if(arp.wtd(x,y)<-1.5){
        arp.transmissivity(x,y) = arp.fdepth(x,y) * arp.ksat(x,y) * \
        std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y));
      }
      // If wtd is greater than 0, max out rate of groundwater movement
      // as though wtd were 0. The surface water will get to move in
      // FillSpillMerge.
      else if(arp.wtd(x,y) > 0){
        arp.transmissivity(x,y) = arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y));
      }
      //Equation S4 from the Fan paper
      else{
        arp.transmissivity(x,y) = arp.ksat(x,y) * (arp.wtd(x,y) + 1.5 + arp.fdepth(x,y));
      }
    }

    // If the fdepth is zero, there is no water transmission below the surface
    // soil layer.
    // If it is less than zero, it is incorrect -- but no water transmission
    // also seems an okay thing to do in this case.
    else{
      arp.transmissivity(x,y) = 0;
    }
  }


  // Updates water-table depth grid over one time step
  #pragma omp parallel for collapse(2) default(none) shared(arp,params)
  for(int32_t y=1; y<params.ncells_y-1; y++)
  for(int32_t x=1; x<params.ncells_x-1; x++){

    // Skip ocean cells
    if(arp.land_mask(x,y) == 1){
      // Otherwise, update the water-table depth change array
      updateCell( params, arp, x, y );;
    }
  }


  // Once all the changes are known, update the WTD everywhere with the
  // difference array
  #pragma omp parallel for collapse(2) default(none) shared(arp, params)
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