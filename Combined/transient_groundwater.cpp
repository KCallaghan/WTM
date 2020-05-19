#include <chrono>
#include <thread>
#include "transient_groundwater.hpp"

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

double FanDarcyGroundwater::computeTransmissivity(ArrayPack &arp, uint32_t x, uint32_t y){
    using namespace std::this_thread;     // sleep_for, sleep_until
    using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
    //cout << "COMPUTING Transmissivity\n";
    double T;
    if(arp.fdepth(x,y)>0){
        // Equation S6 from the Fan paper
        if(arp.wtd(x,y)<-1.5){
            T = arp.fdepth(x,y) * arp.ksat(x,y) \
                       * std::exp( (arp.wtd(x,y)+1.5)/arp.fdepth(x,y) );
        }
        // If wtd is greater than 0, max out rate of groundwater movement
        // as though wtd were 0. The surface water will get to move in
        // FillSpillMerge.
        else if(arp.wtd(x,y) > 0){
            T = arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y));
        }
        //Equation S4 from the Fan paper
        else{
            T = arp.ksat(x,y) * (arp.wtd(x,y) + 1.5 + arp.fdepth(x,y));
        }
    }
    // If the fdepth is zero, there is no water transmission below the surface
    // soil layer.
    // If it is less than zero, it is incorrect -- but no water transmission
    // also seems an okay thing to do in this case.
    else{
        T = 0;
    }
    //if (T != 0){
    //    cout << "T: ";
    //    cout << T;
    //    cout << "\n";
    //    sleep_for(0.3s);
    //}
    return T;
}

void FanDarcyGroundwater::computeNeighborTransmissivity(ArrayPack &arp, uint32_t x, uint32_t y){
    double transmissivityTargetCell = computeTransmissivity(arp, x, y);
    transmissivityN = ( transmissivityTargetCell
                          + computeTransmissivity(arp, x,  y+1) ) / 2.;
    transmissivityS = ( transmissivityTargetCell
                          + computeTransmissivity(arp, x,  y-1) ) / 2.;
    transmissivityW = ( transmissivityTargetCell
                          + computeTransmissivity(arp, x-1,y  ) ) / 2.;
    transmissivityE = ( transmissivityTargetCell
                          + computeTransmissivity(arp, x+1,y  ) ) / 2.;
}

double FanDarcyGroundwater::computeArrayMax(double *_val[], uint8_t size){
    double maxValue = std::numeric_limits<double>::min();
    for (uint32_t i = 0; i < size; i++){
        if(*_val[i] > maxValue){
            maxValue = *_val[i];
        }
    }
    return maxValue;
}

float FanDarcyGroundwater::computeArrayMin(float *_val[], uint8_t size){
    float minValue = std::numeric_limits<double>::max();
    for (uint32_t i = 0; i < size; i++){
        if(*_val[i] < minValue){
            minValue = *_val[i];
        }
    }
    return minValue;
}

double FanDarcyGroundwater::computeMaxStableTimeStep(Parameters &params, ArrayPack &arp, uint32_t x, uint32_t y){
    // Transmissivity is an effective diffusivity
    // Use the highest transmissivity for this time-step calculation
    double dt_max_diffusion_basic;
    double dt_max_diffusion_withPorosity;
    double *Tarray[4] = { &transmissivityN, &transmissivityS,
                             &transmissivityW, &transmissivityE };
    double Dmax = computeArrayMax( Tarray, 4 ); // Max diffusivity
    dt_max_diffusion_basic = ( pow(arp.cellsize_e_w_metres[y], 2)
                                  * pow(params.cellsize_n_s_metres, 2) ) \
                                  / ( 2 * Dmax \
                                    * ( pow(arp.cellsize_e_w_metres[y], 2) \
                                      + pow(params.cellsize_n_s_metres, 2) ) );
    // Now let's add in the porosity differences
    // Porosity differences will AMPLIFY the WTD changes.
    // Let us choose a time step that is also based on the LOWEST porosity
    // (highest WTD change per water volume transfer)
    float *PhiArray[5] = { &arp.porosity(x, y),
                             &arp.porosity(x+1, y), &arp.porosity(x-1, y),
                             &arp.porosity(x, y+1), &arp.porosity(x, y-1) };
    float PhiMin = computeArrayMin( PhiArray, 5 ); // Minimum porosity
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

    //cout << "\n";
    //cout << "DMAX: ";
    //cout << dt_max_diffusion_withPorosity;
    //cout << "\n";
    //cout << "MAX TIME STEP: ";
    //cout << dt_max_diffusion_withPorosity;
    //cout << "\n";
    //cout << "\n";

    return dt_max_diffusion_withPorosity/2.;
}


double FanDarcyGroundwater::computeNewWTD(const float giving_cell_change, const float giving_wtd,const float receiving_wtd, \
	const int x_giving, const int y_giving, const int x_receiving, const int y_receiving, const ArrayPack &arp){

  double receiving_cell_change = 0.0;
  double volume_change = 0.0;
  //we have the change in the giving cell. The giving cell is always losing water. 
  //so, we subtract this value from the giving cell when adjusting wtd later. 

  //first, we check to see if the starting wtd in the giving cell was above the surface. 
  if(giving_wtd > 0){
    //If it stays above 0 once the change has occurred, then no need to worry about porosity. 
    volume_change = giving_cell_change * arp.cell_area[y_giving]; //The volume change is just the height change multiplied by the cell's area. 
    if(giving_wtd - giving_cell_change < 0){  //the water table drops below the surface during this iteration, so we need to consider porosity for part of the water. 
      volume_change = giving_wtd * arp.cell_area[y_giving]; //this is the portion of the water that is above the land surface. 
      volume_change -= arp.cell_area[y_giving] * arp.porosity(x_giving,y_giving) * arp.fdepth(x_giving,y_giving) * \
                      (exp((giving_wtd - giving_cell_change) / arp.fdepth(x_giving,y_giving)) - 1);  //the portion that is below tht land surface. 
        //-= because this comes out as a negative number. 
    }
  }
  else{  // the water table is below the surface to start with, therefore it is below the surface the whole time and we need porosity for all the change.
    volume_change = -arp.cell_area[y_giving] * arp.porosity(x_giving,y_giving) * arp.fdepth(x_giving,y_giving) * \
                      (exp((giving_wtd - giving_cell_change) / arp.fdepth(x_giving,y_giving)) - \
                        exp(giving_wtd / arp.fdepth(x_giving,y_giving)));
  }


  //so now we have the volume change as a positive value from the giving cell, whether it was all above ground, all below ground, or a combination. 
  //Next, we need to use that to calculate the height change in the receiving cell.   

  if(receiving_wtd > 0){  //the receiving cell gains water, so if the starting wtd is above 0, all the change is above the surface. 
    receiving_cell_change = volume_change / arp.cell_area[y_receiving];
  }
  else{  //either it is all below the surface, or a combination. 
 //   //we don't know yet what the height of the change will be, so we start off assuming that it will all be below the surface. 
    receiving_cell_change = arp.fdepth(x_receiving,y_receiving) * log(exp(receiving_wtd / arp.fdepth(x_receiving,y_receiving)) \
          + volume_change / ( arp.cell_area[y_receiving] * arp.porosity(x_receiving,y_receiving) * arp.fdepth(x_receiving,y_receiving)) ) - receiving_wtd;
    if(receiving_wtd +  receiving_cell_change > 0){  //it has changed from GW to SW, so we need to adjust the receiving cell change appropriately. 
      //we want to calculate how much of the water is used up in the ground, i.e. the portion between the starting wtd and 0. 
      double GW_portion = -arp.cell_area[y_receiving] * arp.porosity(x_receiving,y_receiving) * arp.fdepth(x_receiving,y_receiving) * \
                       (exp(receiving_wtd / arp.fdepth(x_receiving,y_receiving)) - 1);
                       //this is the volume of water used up in filling in the ground. 
                       //volume_change - GW_portion is left over to fill surface water. 
      receiving_cell_change = ((volume_change - GW_portion) / arp.cell_area[y_receiving]) - receiving_wtd;
    }

  }
  
  return receiving_cell_change;
}



void FanDarcyGroundwater::computeWTDchangeAtCell(Parameters &params, ArrayPack &arp, int32_t x, int32_t y,
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

    mycell_change = 0.0;

    // Then, compute the discharges
    double QN = transmissivityN * (headN - headCenter) \
                    / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y];
    double QS = transmissivityS * (headS - headCenter) \
                    / params.cellsize_n_s_metres \
                    * arp.cellsize_e_w_metres[y];
    double QE = transmissivityE * (headE - headCenter) \
                    / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres;
    double QW = transmissivityW * (headW - headCenter) \
                    / arp.cellsize_e_w_metres[y] \
                    * params.cellsize_n_s_metres;

    // Move the water, but only in the "internal" variables to handle internal
    // time stepping to maintain stability.
    // dH = sum(discharges) times time step, divided by cell area,
    //      divided by porosity.

  

    wtd_change_N = QN * dt / ( arp.cell_area[y+1]);
    wtd_change_S = QS * dt / ( arp.cell_area[y-1]);
    wtd_change_W = QW * dt / ( arp.cell_area[y] );
    wtd_change_E = QE * dt / ( arp.cell_area[y] );
    


    //Using the wtd_changes from above, we need to calculate how much change will occur in the target cell, accounting for porosity. 

    if(wtd_change_N > 1e-5){  //the current cell is receiving water from the North, so (x,y+1) is the giving cell. 
      //target cell is the receiving cell. 
      mycell_change  += computeNewWTD(wtd_change_N, wtdN, wtdCenter, x, y+1, x, y, arp);
      wtdN +=  - wtd_change_N;
    }
    else if (wtd_change_N < -1e-5){  //the current cell is giving water to the North. The North is the receiving cell. 
      wtdN +=  computeNewWTD(-wtd_change_N, wtdCenter, wtdN, x, y, x, y+1, arp);
      mycell_change  += wtd_change_N;
    }
    else{
      wtd_change_N = 0;
      mycell_change += 0;
      wtdN += 0;
    }

    if(wtd_change_S > 1e-5){
      mycell_change  += computeNewWTD(wtd_change_S, wtdS, wtdCenter, x, y-1, x, y, arp);
      wtdS +=  - wtd_change_S;
    }
    else if (wtd_change_S < -1e-5){
      wtdS +=  computeNewWTD(-wtd_change_S, wtdCenter, wtdS, x, y, x, y-1, arp);
      mycell_change  += wtd_change_S;
    }
    else{
      wtd_change_S = 0;
      mycell_change += 0;
      wtdS += 0;
    }


    if(wtd_change_E > 1e-5){
      mycell_change  += computeNewWTD(wtd_change_E, wtdE, wtdCenter, x+1, y, x, y, arp);
      wtdE +=  - wtd_change_E;
    }
    else if (wtd_change_E < -1e-5){
      wtdE +=  computeNewWTD(-wtd_change_E, wtdCenter, wtdE, x, y, x+1, y, arp);
      mycell_change  += wtd_change_E;
    }
    else{
      wtd_change_E = 0;
      mycell_change += 0;
      wtdE += 0;
    }


    if(wtd_change_W > 1e-5){
      mycell_change  += computeNewWTD(wtd_change_W, wtdW, wtdCenter, x-1, y, x, y, arp);
      wtdW +=  - wtd_change_W;
    }
    else if (wtd_change_W < -1e-5){
      wtdW += computeNewWTD(-wtd_change_W, wtdCenter, wtdW, x, y, x-1, y, arp);
      mycell_change  += wtd_change_W;
    }
    else{
      wtd_change_W = 0;
      mycell_change += 0;
      wtdW += 0;
    }

    //now we have the height changes that will take place in the target cell and each of the four neighbours. 
wtdCenter += mycell_change;

//          if(x==447&&y==715)
  //          std::cout<<"N "<<wtd_change_N<<" S "<<wtd_change_S<<" E "<<wtd_change_E<<" W "<<wtd_change_W<<std::endl;
}

void FanDarcyGroundwater::updateCell(Parameters &params, ArrayPack &arp, uint32_t x, uint32_t y){

    using namespace std::this_thread;     // sleep_for, sleep_until
    using namespace std::chrono_literals; // ns, us, ms, s, h, etc.

    // Runs functions to compute time steps and update WTD for the center cell
    // and its neighbours until the outer time step has been completed

    // Initialize variables for dynamic time stepping
    double time_remaining = params.deltat;
    double dt_inner;

    // Initial water-table depths, prior to updating
    double wtdCenter_initial = arp.wtd(x,y);
    wtdCenter = arp.wtd(x,y);


    double temp = params.deltat;
    //cout << "\n";
    //cout << "\n";
    //cout << wtdCenter_initial;
    //cout << "\n";
    wtdN      = arp.wtd(x,y+1);
    wtdS      = arp.wtd(x,y-1);
    wtdE      = arp.wtd(x+1,y);
    wtdW      = arp.wtd(x-1,y);
    //cout << wtdW;
    //cout << "\n";
    //cout << "\n";

    // Update water-table depths using dynamic time stepping
    while (time_remaining > 0){
        //cout << wtdCenter;
        //cout << " ";
        //cout << time_remaining;
        //cout << "\n";
        computeNeighborTransmissivity(arp, x, y);  //currently transmissivity is based on wtd, which does not change during the while loop.
        //should we change it during the while loop? If not, we can move it out. 


        for(int32_t y=1; y<params.ncells_y-1; y++){
        for(int32_t x=1; x<params.ncells_x-1; x++){
 
        double max_stable_time_step = computeMaxStableTimeStep(params, arp, x, y);

        if(max_stable_time_step<temp)
        	temp = max_stable_time_step;

    }
}
        // Choose the inner-loop time step
        if(time_remaining <= temp){
            dt_inner = time_remaining;
        }
        else{
            dt_inner = temp;
        }
        //cout << "\n";
        //cout << "\n";
        //cout << "!dt_inner\n";
        //cout << dt_inner;
        //cout << "\n";
        //cout << "!dt_inner\n";
        //cout << "\n";
        //cout << "\n";
        computeWTDchangeAtCell(params, arp, x, y, dt_inner);
        time_remaining -= dt_inner;
        //cout << wtdCenter;
        //cout << "\n";
        //cout << arp.cell_area[y];
        //cout << "\n";
        //cout << arp.porosity(x,y);
        //cout << "\n";
        //cout << params.cellsize_n_s_metres;
        //cout << "\n";
        //cout << arp.cellsize_e_w_metres[y];
        //cout << "\n";
        //sleep_for(0.3s);
 //       if(x==447&&y==715)
   //     	std::cout<<"dt inner was "<<dt_inner<<" change so far was "<<wtdCenter - wtdCenter_initial<<std::endl;
    }
    // When exiting loop, the wtdCenter variable holds the final
    // water-table depth
    // This subtraction is unnecessary; could just give the new total instead
    arp.wtd_change_total(x,y) = wtdCenter - wtdCenter_initial;
    //cout << "wtd_change_total: ";
    //cout << arp.wtd_change_total(x,y);
    //cout << " m\n";
}

/////////////////
// CONSTRUCTOR //
/////////////////

FanDarcyGroundwater::FanDarcyGroundwater(){
}

//FanDarcyGroundwater::FanDarcyGroundwater(Parameters _params, ArrayPack _arp){
//    arp = _arp;
//    params = _params;
//}

//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

//void FanDarcyGroundwater::set_arp(ArrayPack _arp){
//    arp = _arp;
//}

//void FanDarcyGroundwater::set_params(Parameters _params){
//    params = _params;
//}

void FanDarcyGroundwater::initialize(){

}

void FanDarcyGroundwater::update(Parameters &params, ArrayPack &arp){
    //cout << "Output sentence\n"; // prints Output sentence on screen
    //cout << params.ncells_x; // prints Output sentence on screen
    //cout << "\n";
    //cout << params.ncells_y; // prints Output sentence on screen
    //cout << "\n";
    // Updates water-table depth grid over one time step
    for(int32_t y=1; y<params.ncells_y-1; y++){
        for(int32_t x=1; x<params.ncells_x-1; x++){
            // Skip ocean cells
            if(arp.land_mask(x,y) == 0){
                continue;
            }
            // Otherwise, update the water-table depth change array
            updateCell( params, arp, x, y );
            //cout << "Test";
            //cout << arp.wtd_change_total(x,y);
            //cout << "\n";
        }
    }
    // Once all the changes are known, update the WTD everywhere with the
    // difference array
   // cout << arp.wtd(447,715);
   // cout << "-->";
    for(int32_t y=1; y<params.ncells_y-1; y++){
        for(int32_t x=1; x<params.ncells_x-1; x++){
            // Skip ocean cells
            if(arp.land_mask(x,y) == 0){
                continue;
            }
            // Update the whole wtd array at once.
            // This is the new water table after groundwater has moved
            // for delta_t seconds.
            //cout << arp.wtd(x,y);
            //cout << "-->";
            arp.wtd(x,y) += arp.wtd_change_total(x,y);
            //cout << arp.wtd(x,y);
            //cout << "\n";
        }
    }
  //  cout << arp.wtd(447,715);
  //  cout << "\n";
    //SaveAsNetCDF(arp.wtd, "wtdCheck.nc", "WTD");
}

// This can be populated if we intend to run this module on its own.
// Otherwise, will not be called
void FanDarcyGroundwater::run(){

}

// This can include functions to clean up / clear memory, if needed
void FanDarcyGroundwater::finalize(){

}
