#include "transient_groundwater.hpp"

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

float64_t FanDarcyGroundwater::computeTransmissivity(int x, int y){
    if(arp.fdepth(x,y)>0){
        // Equation S6 from the Fan paper
        if(arp.wtd(x,y)<-1.5){
            return arp.fdepth(x,y) * arp.ksat(x,y) \
                       * std::exp( (arp.wtd(x,y)+1.5)/arp.fdepth(x,y) );
        }
        // If wtd is greater than 0, max out rate of groundwater movement
        // as though wtd were 0. The surface water will get to move in
        // FillSpillMerge.
        else if(arp.wtd(x,y) > 0){
            return arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y));
        }
        //Equation S4 from the Fan paper
        else{
            return arp.ksat(x,y) * (arp.wtd(x,y) + 1.5 + arp.fdepth(x,y));
        }
    }
    // If the fdepth is zero, there is no water transmission below the surface
    // soil layer.
    // If it is less than zero, it is incorrect -- but no water transmission
    // also seems an okay thing to do in this case.
    else{
        return 0;
    }
}

void FanDarcyGroundwater::computeNeighborTransmissivity(int32_t x, int32_t y){
    float64_t transmissivityTargetCell = computeTransmissivity(x, y);
    transmissivityN = ( transmissivityTargetCell + kcell(x,  y+1) ) / 2.;
    transmissivityS = ( transmissivityTargetCell + kcell(x,  y-1) ) / 2.;
    transmissivityW = ( transmissivityTargetCell + kcell(x-1,y  ) ) / 2.;
    transmissivityE = ( transmissivityTargetCell + kcell(x+1,y  ) ) / 2.;
}

float64_t FanDarcyGroundwater::computeArrayMax(float64_t T, uint32_t size){
    float64_t maxValue = std::numeric_limits<float64_t>::min();
    for (i = 0; i < size; i++){
        if(T[i] > maxValue){
            maxValue = T[i];
        }
    }
    return maxValue;
}

float64_t FanDarcyGroundwater::computeArrayMin(float64_t T, uint32_t size){
    float64_t minValue = std::numeric_limits<float64_t>::max();
    for (i = 0; i < size; i++){
        if(T[i] < minValue){
            minValue = T[i];
        }
    }
    return minValue;
}

float64_t FanDarcyGroundwater::computeMaxStableTimeStep(int32_t x, int32_t y){
    // Transmissivity is an effective diffusivity
    // Use the highest transmissivity for this time-step calculation
    float64_t dt_max_diffusion_basic;
    float64_t dt_max_diffusion_withPorosity;
    float64_t Tarray[4] = { &transmissivityN, &transmissivityS,
                            &transmissivityW, &transmissivityE };
    float64_t Dmax = computeArrayMax( Tarray, 4 ); // Max diffusivity
    dt_max_diffusion_basic = ( arp.cellsize_e_w_metres[y]**2
                                  * params.cellsize_n_s_metres**2) \
                                  / ( 2 * Dmax \
                                      * ( arp.cellsize_e_w_metres[y]**2
                                          + arp.cellsize_e_w_metres[y]**2 ) );
    // Now let's add in the porosity differences
    // Porosity differences will AMPLIFY the WTD changes.
    // Let us choose a time step that is also based on the LOWEST porosity
    // (highest WTD change per water volume transfer)
    float64_t PhiArray[5] = { &arp.porosity(x, y),
                              &arp.porosity(x+1, y), &arp.porosity(x-1, y),
                              &arp.porosity(x, y+1), &arp.porosity(x, y-1) };
    float64_t PhiMin = computeArrayMin( Tarray, 4 ); // Minimum porosity
    // Porosity is a linear amplifier of WTD change, and it amplifies change
    // in both the giving and receiving cells.
    // Amplification goes as 1 / phi.
    // We need to consider this going up and down, so 1 / (2*phi)
    dt_max_diffusion_withPorosity = dt_max_diffusion_basic / (2 * PhiMin);
    // In order to avoid operating at the very maximum time step possible,
    // we apply a factor of safety of 2
    return dt_max_diffusion_withPorosity/2.;
}

void FanDarcyGroundwater::computeCellNeighborDischarge(int32_t x, int32_t y){
    // Get the hydraulic conductivity for our cells of interest

}



void receiving_cell_wtd(float giving_cell_change, float giving_wtd,
                          float receiving_wtd, int x_giving, int y_giving,
                          int x_receiving, int y_receiving){
}




double get_change(const int x, const int y, const double time_remaining){

  // Declare variables updated in the loop
  double wtd_change_N;
  double wtd_change_S;
  double wtd_change_E;
  double wtd_change_W;

  // Elevation head - topography plus the water table depth (negative if
  // water table is below earth surface)
  const auto my_head = arp.topo(x,y) + params.wtdCenter;
  // heads for each of my neighbour cells
  const auto headN   = arp.topo(x,y+1) + params.wtdN;
  const auto headS   = arp.topo(x,y-1) + params.wtdS;
  const auto headW   = arp.topo(x-1,y) + params.wtdW;
  const auto headE   = arp.topo(x+1,y) + params.wtdE;


  float mycell_change_N = 0.0;
  float mycell_change_S = 0.0;
  float mycell_change_E = 0.0;
  float mycell_change_W = 0.0;
  float mycell_change   = 0.0;

  float change_in_N_cell = 0.0;
  float change_in_S_cell = 0.0;
  float change_in_E_cell = 0.0;
  float change_in_W_cell = 0.0;

  bool stable = false;
  double time_step = time_remaining;

  //if( x== 1366&& y ==794)
  // std::cout<<"************************************"<<std::endl;

  while(!stable){

    //if( x== 1366&& y ==794)
    // std::cout<<"params S "<<params.S<<" wtd change S "<<wtd_change_S<<" change "<<change_in_S_cell<<" time "<<time_step<<" remaining "<<time_remaining<<std::endl;

    //   std::cout<<"time step is "<<time_step<<" and time remaining is "<<time_remaining<<std::endl;
    // Change in water-table depth.
    // (1) Discharge across cell boundaries
    // Average hydraulic conductivity of the two cells *
    // head difference between the two / distance (i.e., dH/dx_i) *
    // width of cell across which the water is discharged *
    // time step
    // (2) Divide by the area of the given cell: maps water volume
    // increases/decreases to change in head
    //Q_N = kN * (headN - my_head)
    //Q_S =
    //Q_W =
    //Q_E =


    wtd_change_N = kN * (headN - my_head) / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y] * time_step \
                    / arp.cell_area[y];
    wtd_change_S = kS * (headS - my_head) / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y] * time_step \
                    / arp.cell_area[y];
    wtd_change_E = kE * (headE - my_head) / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres * time_step \
                    / arp.cell_area[y];
    wtd_change_W = kW * (headW - my_head) / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres * time_step \
                    / arp.cell_area[y];


    //Using the wtd_changes from above, we need to calculate how much change will occur in the target cell, accounting for porosity.

    if(wtd_change_N > 1e-5){  //the current cell is receiving water from the North, so (x,y+1) is the giving cell.
      //target cell is the receiving cell.
      change_in_N_cell =  - wtd_change_N;
      mycell_change_N  = receiving_cell_wtd(wtd_change_N, params.N, params.me, x, y+1, x, y, arp);
    }
    else if (wtd_change_N < -1e-5){  //the current cell is giving water to the North. The North is the receiving cell.
      change_in_N_cell =  receiving_cell_wtd(-wtd_change_N, params.me, params.N, x, y, x, y+1, arp);
      mycell_change_N  = wtd_change_N;
    }
    else{
      wtd_change_N = 0;
      mycell_change_N = 0;
      change_in_N_cell = 0;
    }

    if(wtd_change_S > 1e-5){
      change_in_S_cell =  - wtd_change_S;
      mycell_change_S  = receiving_cell_wtd(wtd_change_S, params.S, params.me, x, y-1, x, y, arp);
    }
    else if (wtd_change_S < -1e-5){
      change_in_S_cell =  receiving_cell_wtd(-wtd_change_S, params.me, params.S, x, y, x, y-1, arp);
      mycell_change_S  = wtd_change_S;
    }
    else{
      wtd_change_S = 0;
      mycell_change_S = 0;
      change_in_S_cell = 0;
    }

    if(wtd_change_E > 1e-5){
      change_in_E_cell =  - wtd_change_E;
      mycell_change_E  = receiving_cell_wtd(wtd_change_E, params.E, params.me, x+1, y, x, y, arp);
    }
    else if (wtd_change_E < -1e-5){
      change_in_E_cell =  receiving_cell_wtd(-wtd_change_E, params.me, params.E, x, y, x+1, y, arp);
      mycell_change_E  = wtd_change_E;
    }
    else{
      wtd_change_E = 0;
      mycell_change_E = 0;
      change_in_E_cell = 0;
    }

    if(wtd_change_W > 1e-5){
      change_in_W_cell =  - wtd_change_W;
      mycell_change_W  = receiving_cell_wtd(wtd_change_W, params.W, params.me, x-1, y, x, y, arp);
    }
    else if (wtd_change_W < -1e-5){
      change_in_W_cell = receiving_cell_wtd(-wtd_change_W, params.me, params.W, x, y, x-1, y, arp);
      mycell_change_W  = wtd_change_W;
    }
    else{
      wtd_change_W = 0;
      mycell_change_W = 0;
      change_in_W_cell = 0;
    }

    //now we have the height changes that will take place in the target cell and each of the four neighbours.

    //Total change in wtd for our target cell in this iteration

    mycell_change =  mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W;

    if( ((headN > my_head) && (headS > my_head) && (change_in_N_cell+arp.topo(x,y+1)+params.N < mycell_change+arp.topo(x,y)+params.me) && (change_in_S_cell+arp.topo(x,y-1)+params.S < mycell_change+arp.topo(x,y)+params.me) && fabs(wtd_change_N)> 1e-6 && fabs(wtd_change_S) > 1e-6 ) ||  \
        ((headN < my_head) && (headS < my_head) && (change_in_N_cell+arp.topo(x,y+1)+params.N > mycell_change+arp.topo(x,y)+params.me) && (change_in_S_cell+arp.topo(x,y-1)+params.S > mycell_change+arp.topo(x,y)+params.me) && fabs(wtd_change_N)> 1e-6 && fabs(wtd_change_S) > 1e-6 ) ||  \
        ((headE > my_head) && (headW > my_head) && (change_in_E_cell+arp.topo(x+1,y)+params.E < mycell_change+arp.topo(x,y)+params.me) && (change_in_W_cell+arp.topo(x-1,y)+params.W < mycell_change+arp.topo(x,y)+params.me) && fabs(wtd_change_E)> 1e-6 && fabs(wtd_change_W) > 1e-6 ) ||  \
        ((headE < my_head) && (headW < my_head) && (change_in_E_cell+arp.topo(x+1,y)+params.E > mycell_change+arp.topo(x,y)+params.me) && (change_in_W_cell+arp.topo(x-1,y)+params.W > mycell_change+arp.topo(x,y)+params.me) && fabs(wtd_change_E)> 1e-6 && fabs(wtd_change_W) > 1e-6 )  ){
      //there is an instability.
      time_step = time_step/2.;
    }
    else if( (((headN - my_head)*((headN + change_in_N_cell) - (my_head + mycell_change_N)) < 0) && wtd_change_N > 1e-6)  || \
             (((headS - my_head)*((headS + change_in_S_cell) - (my_head + mycell_change_S)) < 0) && wtd_change_S > 1e-6) ||
             (((headE - my_head)*((headE + change_in_E_cell) - (my_head + mycell_change_E)) < 0) && wtd_change_E > 1e-6) ||
             (((headW - my_head)*((headW + change_in_W_cell) - (my_head + mycell_change_W)) < 0) && wtd_change_W > 1e-6) ){
      // The change between any 2 cells can't be greater than the difference
      // between those two cells.
      time_step = time_step/2.;
    }
    else{
      // Otherwise, a stable time step has been found
      arp.wtd_change_total(x,y) += ( mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W );
      params.me += mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W;
      params.N  += change_in_N_cell;
      params.S  += change_in_S_cell;
      params.W  += change_in_W_cell;
      params.E  += change_in_E_cell;

      stable = true;

      if(fabs(arp.wtd_change_total(x,y)) > 100000){
        std::cout<<"stable "<<arp.wtd(x,y)<<" change "<<arp.wtd_change_total(x,y)<<" x "<<x<<" y "<<y<<std::endl;

      std::cout<<"N "<<mycell_change_N<<" S "<<mycell_change_S<<" E "<<mycell_change_E<<" W "<<mycell_change_W<<std::endl;
      std::cout<<"wtd change N "<<wtd_change_N<<" S "<<wtd_change_S<<" E "<<wtd_change_E<<" W "<<wtd_change_W<<std::endl;
      std::cout<<"wtd W "<<arp.wtd(x-1,y)<<" wtd me "<<arp.wtd(x,y)<<std::endl;
      std::cout<<"W was the problem "<<headW<<" mine "<<my_head<<" k "<<kW<<" time step "<<time_step<<std::endl;

      //wtd_change_S = kS * (headS - my_head) / params.cellsize_n_s_metres \
                          * arp.cellsize_e_w_metres[y] * time_step \
                          / arp.cell_area[y];
      }
    }
  }
  return time_step;
}




void groundwater(){
  /**
  @param params   Global paramaters - we use the texfilename, run type,
                  number of cells in the x and y directions (ncells_x
                  and ncells_y), delta_t (number of seconds in a time step),
                  and cellsize_n_s_metres (size of a cell in the north-south
                  direction)

  @param arp      Global arrays - we access land_mask, topo, wtd,
                  wtd_change_total, cellsize_e_w_metres, and cell_area.
                  land_mask is a binary representation of where land is vs
                                        where ocean is.
                  topo is the input topography, i.e. land elevation above
                                        sea level in each cell.
                  wtd is the water table depth.
                  wtd_change_total is the amount by which wtd will change
                                        during this time step as a result of
                                        groundwater movement.
                  cellsize_e_w_metres is the distance across a cell in the
                                        east-west direction.
                  cell_area is the area of the cell, needed because different
                                        cells have different areas and
                                        therefore an accommodate different
                                        water volumes.

  @return  An updated wtd that represents the status of the water table after
           groundwater has been able to flow for the amount of time represented
           by delta_t.
  **/


  // Declare status variables and set initial values to 0
  double total_changes = 0.;
  float max_total      = 0.;
  float min_total      = 0.;
  float max_change     = 0.;


  // Set up log file
  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);

  textfile<<"Groundwater"<<std::endl;


  ////////////////////////////
  // COMPUTE CHANGES IN WTD //
  ////////////////////////////

  // Cycle through the entire array, calculating how much the water-table
  // changes in each cell per iteration.
  // We do this instead of using a staggered grid to approx. double CPU time
  // in exchange for using less memory.
  for(int y=1; y<params.ncells_y-1; y++){
    for(int x=1; x<params.ncells_x-1; x++){
      //skip ocean cells

      if(arp.land_mask(x,y) == 0)
        continue;

      double time_step = 0.0;
      arp.wtd_change_total(x,y) = 0.0;
      double time_remaining = params.deltat;

      params.N  = arp.wtd(x,y+1);
      params.S  = arp.wtd(x,y-1);
      params.E  = arp.wtd(x+1,y);
      params.W  = arp.wtd(x-1,y);
      params.me = arp.wtd(x,y);

      while(time_remaining > 1e-4){
        time_step = get_change(x, y,time_remaining,params,arp);
        time_remaining -= time_step;
      }
    }
  }


  ////////////////
  // UPDATE WTD //
  ////////////////

  for(int y=1;y<params.ncells_y-1;y++){
    for(int x=1;x<params.ncells_x-1; x++){

      if(arp.land_mask(x,y) == 0)
        continue;

      // Update the whole wtd array at once.
      // This is the new water table after groundwater has moved
      // for delta_t seconds.
      arp.wtd(x,y) += arp.wtd_change_total(x,y);
    }
  }

  // Write status to text file
  textfile << "total GW changes were " << total_changes << std::endl;
  textfile << "max wtd was " << max_total << " and min wtd was " \
           << min_total << std::endl;
  textfile << "max GW change was " << max_change << std::endl;
  textfile.close();
}
