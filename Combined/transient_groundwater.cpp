#include "transient_groundwater.hpp"

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

double FanDarcyGroundwater::computeTransmissivity(uint32_t x, uint32_t y){
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

void FanDarcyGroundwater::computeNeighborTransmissivity(uint32_t x, uint32_t y){
    double transmissivityTargetCell = computeTransmissivity(x, y);
    transmissivityN = ( transmissivityTargetCell
                          + computeTransmissivity(x,  y+1) ) / 2.;
    transmissivityS = ( transmissivityTargetCell
                          + computeTransmissivity(x,  y-1) ) / 2.;
    transmissivityW = ( transmissivityTargetCell
                          + computeTransmissivity(x-1,y  ) ) / 2.;
    transmissivityE = ( transmissivityTargetCell
                          + computeTransmissivity(x+1,y  ) ) / 2.;
}

double FanDarcyGroundwater::computeArrayMax(double T[], uint8_t size){
    double maxValue = std::numeric_limits<double>::min();
    for (uint32_t i = 0; i < size; i++){
        if(T[i] > maxValue){
            maxValue = T[i];
        }
    }
    return maxValue;
}

double FanDarcyGroundwater::computeArrayMin(double T[], uint8_t size){
    double minValue = std::numeric_limits<double>::max();
    for (uint32_t i = 0; i < size; i++){
        if(T[i] < minValue){
            minValue = T[i];
        }
    }
    return minValue;
}

double FanDarcyGroundwater::computeMaxStableTimeStep(uint32_t x, uint32_t y){
    // Transmissivity is an effective diffusivity
    // Use the highest transmissivity for this time-step calculation
    double dt_max_diffusion_basic;
    double dt_max_diffusion_withPorosity;
    double Tarray[4] = { &transmissivityN, &transmissivityS,
                            &transmissivityW, &transmissivityE };
    double Dmax = computeArrayMax( Tarray, 4 ); // Max diffusivity
    dt_max_diffusion_basic = ( arp.cellsize_e_w_metres[y]**2
                                  * params.cellsize_n_s_metres**2) \
                                  / ( 2 * Dmax \
                                      * ( arp.cellsize_e_w_metres[y]**2
                                          + arp.cellsize_e_w_metres[y]**2 ) );
    // Now let's add in the porosity differences
    // Porosity differences will AMPLIFY the WTD changes.
    // Let us choose a time step that is also based on the LOWEST porosity
    // (highest WTD change per water volume transfer)
    double PhiArray[5] = { &arp.porosity(x, y),
                              &arp.porosity(x+1, y), &arp.porosity(x-1, y),
                              &arp.porosity(x, y+1), &arp.porosity(x, y-1) };
    double PhiMin = computeArrayMin( Tarray, 4 ); // Minimum porosity
    // Porosity is a linear amplifier of WTD change, and it amplifies change
    // in both the giving and receiving cells.
    // Amplification goes as 1 / phi.
    // We need to consider this going up and down, so 1 / (2*phi)
    dt_max_diffusion_withPorosity = dt_max_diffusion_basic / (2 * PhiMin);
    // In order to avoid operating at the very maximum time step possible,
    // we apply a factor of safety of 2
    return dt_max_diffusion_withPorosity/2.;
}

void FanDarcyGroundwater::computeWTDchangeAtCell(int32_t x, int32_t y,
                                                 double dt){
    // Update WTD change

    // We do this instead of using a staggered grid to approx. double CPU time
    // in exchange for using less memory.

    // First, compute elevation head at center cell and all neighbours
    // This equals topography plus the water table depth
    // (positive upwards; negative if water table is below Earth surface)
    double headCenter = arp.topo(x,y) + wtdCenter;
    double headN      = arp.topo(x,y+1) + wtdN;
    double headS      = arp.topo(x,y-1) + wtdS;
    double headW      = arp.topo(x-1,y) + wtdW;
    double headE      = arp.topo(x+1,y) + wtdE;

    // Then, compute the discharges
    double QN = transmissivityN * (headN - my_head) \
                    / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y];
    double QS = transmissivityS * (headS - my_head) \
                    / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y];
    double QE = transmissivityE * (headS - my_head) \
                    / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres;
    double QW = transmissivityW * (headW - my_head) \
                    / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres;

    // Move the water, but only in the "internal" variables to handle internal
    // time stepping to maintain stability.
    // dH = sum(discharges) times time step, divided by cell area,
    //      divided by porosity.
    wtdCenter = ( QN + QS + QE + QW ) * dt \
                           / ( arp.cell_area[y] * arp.porosity(x,y) )
    wtdN -= QN * dt / ( arp.cell_area[y+1] * arp.porosity(x,y+1) )
    wtdS -= QS * dt / ( arp.cell_area[y-1] * arp.porosity(x,y-1) )
    wtdW -= QW * dt / ( arp.cell_area[y-1] * arp.porosity(x,y-1) )
    wtdE -= QE * dt / ( arp.cell_area[y+1] * arp.porosity(x,y+1) )
}

void FanDarcyGroundwater::updateCell(uint32_t x, uint32_t y){
    // Runs functions to compute time steps and update WTD for the center cell
    // and its neighbours until the outer time step has been completed

    // Initialize variables for dynamic time stepping
    double time_remaining = params.deltat;
    double dt_inner;

    // Initial water-table depths, prior to updating
    wtdCenter = arp.wtd(x,y);
    wtdN      = arp.wtd(x,y+1);
    wtdS      = arp.wtd(x,y-1);
    wtdE      = arp.wtd(x+1,y);
    wtdW      = arp.wtd(x-1,y);

    // Update water-table depths using dynamic time stepping
    while (time_remaining > 0){
        double max_stable_time_step = computeMaxStableTimeStep(x,y);
        // Choose the inner-loop time step
        if(time_remaining >= max_stable_time_step){
            dt_inner = time_remaining;
        }
        else{
            dt_inner = max_stable_time_step;
        }
        computeWTDchangeAtCell(x, y, dt_inner);
        time_remaining -= dt_inner;
    }
    // When exiting loop, the wtdCenter variable holds the final
    // water-table depth
    arp.wtd_change_total(x,y) = wtdCenter;
}

//!WHERE IS TOTAL_CHANGES UPDATED???????????????????????????????????????????????
void FanDarcyGroundwater::logToFile(){
  // Set up log file
  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);
  textfile<<"Groundwater"<<std::endl;
  textfile << "total GW changes were " << total_changes << std::endl;
  textfile << "max wtd was " << max_total << " and min wtd was " \
           << min_total << std::endl;
  textfile << "max GW change was " << max_change << std::endl;
  textfile.close();
}

/////////////////
// CONSTRUCTOR //
/////////////////

FanDarcyGroundwater::FanDarcyGroundwater(){
}

//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void set_arp(ArrayPack &_arp){
    arp = &_arp;
}

void set_params(Parameters &_params){
    params = &_params;
}


void FanDarcyGroundwater::initialize(){

}

void FanDarcyGroundwater::update(bool _log){
    // Updates water-table depth grid over one time step
    for(uint32_t y=1; y<params.ncells_y-1; y++){
        for(uint32_t x=1; x<params.ncells_x-1; x++){
            // Skip ocean cells
            if(arp.land_mask(x,y) == 0);
                continue;
            // Otherwise, update the water-table depth change array
            updateCell( x, y );
        }
    }
    // Once all the changes are known, update the WTD everywhere with the
    // difference array
    for(int y=1;y<params.ncells_y-1;y++){
        for(int x=1;x<params.ncells_x-1; x++){
            // Skip ocean cells
            if(arp.land_mask(x,y) == 0){
                continue;
            }
            // Update the whole wtd array at once.
            // This is the new water table after groundwater has moved
            // for delta_t seconds.
            arp.wtd(x,y) += arp.wtd_change_total(x,y);
        }
    }
    if(_log){
        logToFile();
    }
}

// This can be populated if we intend to run this module on its own.
// Otherwise, will not be called
void FanDarcyGroundwater::run(){

}

// This can include functions to clean up / clear memory, if needed
void FanDarcyGroundwater::finalize(){

}
