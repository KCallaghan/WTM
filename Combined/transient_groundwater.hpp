#include "../common/netcdf.hpp"
#include "ArrayPack.hpp"
#include "parameters.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
using namespace std;

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;

double kcell(const int x, const int y, const ArrayPack &arp){
  /**
  Mini-function that gives the hydraulic conductivity per cell, kcell.
  This changes through time as the water-table depth changes:

  @param x         The x-coordinate of the cell in question

  @param y         The y-coordinate of the cell in question

  @param ArrayPack Global arrays. Here we use: 
          - fdepth: The e-folding depth based on slope and temperature. 
                    This describes the decay of kcell with depth. 
          - wtd:    The water-table depth. We use a different calculation
                    For water tables above vs below 1.5 m below land surface.
          - ksat:   Hydraulic conductivity, based on soil types
  @return  The kcell value for the cell in question. This is the 
           integration of the hydraulic conductivity over flow depth.
  **/

  if(arp.fdepth(x,y)>0){
    // Equation S6 from the Fan paper
    if(arp.wtd(x,y)<-1.5)
      return arp.fdepth(x,y) * arp.ksat(x,y) \
                 * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); 
    // If wtd is greater than 0, max out rate of groundwater movement 
    // as though wtd were 0. The surface water will get to move in 
    // FillSpillMerge.
    else if(arp.wtd(x,y) > 0)
      return arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y));                                        
    else
      return arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));                             //Equation S4 from the Fan paper
  }else
    return 0;
}





double receiving_cell_wtd(const float giving_cell_change, const float giving_wtd,const float receiving_wtd, const int x_giving, const int y_giving, \
                          const int x_receiving, const int y_receiving, const ArrayPack &arp){

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
    //we don't know yet what the height of the change will be, so we start off assuming that it will all be below the surface. 
    receiving_cell_change =  arp.fdepth(x_receiving,y_receiving) * log(exp(receiving_wtd / arp.fdepth(x_receiving,y_receiving)) \
          + volume_change / (arp.cell_area[y_receiving] * arp.porosity(x_receiving,y_receiving) * arp.fdepth(x_receiving,y_receiving)) ) - receiving_wtd;
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



double get_change(const int x, const int y, const double time_remaining, Parameters &params, ArrayPack &arp){
  
  // Declare variables updated in the loop
  double wtd_change_N;
  double wtd_change_S;
  double wtd_change_E;
  double wtd_change_W;


  // Elevation head - topography plus the water table depth (negative if 
  // water table is below earth surface)               
  const auto my_head = arp.topo(x,y) + params.me;
  // heads for each of my neighbour cells
  const auto headN   = arp.topo(x,y+1) + params.N;
  const auto headS   = arp.topo(x,y-1) + params.S;              
  const auto headW   = arp.topo(x-1,y) + params.W;              
  const auto headE   = arp.topo(x+1,y) + params.E;              

  // Get the hydraulic conductivity for our cells of interest
  const auto kN = ( kcell(x,  y,   arp) + kcell(x,  y+1, arp) ) / 2.;
  const auto kS = ( kcell(x,  y,   arp) + kcell(x,  y-1, arp) ) / 2.;                
  const auto kW = ( kcell(x,  y,   arp) + kcell(x-1,y,   arp) ) / 2.;
  const auto kE = ( kcell(x,  y,   arp) + kcell(x+1,y,   arp) ) / 2.;

  float mycell_change_N = 0.0;
  float mycell_change_S = 0.0;
  float mycell_change_E = 0.0;
  float mycell_change_W = 0.0;
  float mycell_change   = 0.0;

  float change_in_N_cell = 0.0;
  float change_in_S_cell = 0.0;
  float change_in_E_cell = 0.0;
  float change_in_W_cell = 0.0;

  int stable = 0;
  double time_step = time_remaining;

  while(stable == 0){
 //   std::cout<<"time step is "<<time_step<<" and time remaining is "<<time_remaining<<std::endl;
    // Change in water-table depth.
    // (1) Discharge across cell boundaries
    // Average hydraulic conductivity of the two cells * 
    // head difference between the two / distance (i.e., dH/dx_i) *
    // width of cell across which the water is discharged *
    // time step
    // (2) Divide by the area of the given cell: maps water volume 
    // increases/decreases to change in head
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

    if(wtd_change_N > 0){  //the current cell is receiving water from the North, so (x,y+1) is the giving cell. 
      //target cell is the receiving cell. 
      change_in_N_cell = arp.topo(x,y+1) + params.N - wtd_change_N;
      mycell_change_N  = receiving_cell_wtd(wtd_change_N, params.N, params.me, x, y+1, x, y, arp);
    }
    else{  //the current cell is giving water to the North. The North is the receiving cell. 
      change_in_N_cell = arp.topo(x,y+1) + params.N + receiving_cell_wtd(-wtd_change_N, params.me, params.N, x, y, x, y+1, arp);
      mycell_change_N  = wtd_change_N;
    }

    if(wtd_change_S > 0){
      change_in_S_cell = arp.topo(x,y-1) + params.S - wtd_change_S;
      mycell_change_S  = receiving_cell_wtd(wtd_change_S, params.S, params.me, x, y-1, x, y, arp);
    }
    else{
      change_in_S_cell = arp.topo(x,y-1) + params.S + receiving_cell_wtd(-wtd_change_S, params.me, params.S, x, y, x, y-1, arp);
      mycell_change_S  = wtd_change_S;
    }

    if(wtd_change_E > 0){
      change_in_E_cell = arp.topo(x+1,y) + params.E - wtd_change_E;
      mycell_change_E  = receiving_cell_wtd(wtd_change_E, params.E, params.me, x+1, y, x, y, arp);
    }
    else{
      change_in_E_cell = arp.topo(x+1,y) + params.E + receiving_cell_wtd(-wtd_change_E, params.me, params.E, x, y, x+1, y, arp);
      mycell_change_E  = wtd_change_E;
    }

    if(wtd_change_W > 0){
      change_in_W_cell = arp.topo(x-1,y) + params.W - wtd_change_W;
      mycell_change_W  = receiving_cell_wtd(wtd_change_W, params.W, params.me, x-1, y, x, y, arp);
    }
    else{
      change_in_W_cell = arp.topo(x-1,y) + params.W + receiving_cell_wtd(-wtd_change_W, params.me, params.W, x, y, x-1, y, arp);
      mycell_change_W  = wtd_change_W;
    }
    //now we have the height changes that will take place in the target cell and each of the four neighbours. 

    //Total change in wtd for our target cell in this iteration


    mycell_change = arp.topo(x,y) + params.me + mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W;

    if( ((headN > my_head) && (headS > my_head) && (change_in_N_cell < mycell_change) && (change_in_S_cell < mycell_change) && fabs(wtd_change_N)> 1e-6 && fabs(wtd_change_S) > 1e-6 ) ||  \
        ((headN < my_head) && (headS < my_head) && (change_in_N_cell > mycell_change) && (change_in_S_cell > mycell_change) && fabs(wtd_change_N)> 1e-6 && fabs(wtd_change_S) > 1e-6 ) ||  \
        ((headE > my_head) && (headW > my_head) && (change_in_E_cell < mycell_change) && (change_in_W_cell < mycell_change) && fabs(wtd_change_E)> 1e-6 && fabs(wtd_change_W) > 1e-6 ) ||  \
        ((headE < my_head) && (headW < my_head) && (change_in_E_cell > mycell_change) && (change_in_W_cell > mycell_change) && fabs(wtd_change_E)> 1e-6 && fabs(wtd_change_W) > 1e-6 )  ){
      //there is an instability. 
      time_step = time_step/2.;
    }
    else if( (((headN - my_head)*(change_in_N_cell - (my_head + mycell_change_N)) < 0) && wtd_change_N > 1e-6)  || \ 
             (((headS - my_head)*(change_in_S_cell - (my_head + mycell_change_S)) < 0) && wtd_change_S > 1e-6) ||
             (((headE - my_head)*(change_in_E_cell - (my_head + mycell_change_E)) < 0) && wtd_change_E > 1e-6) ||
             (((headW - my_head)*(change_in_W_cell - (my_head + mycell_change_W)) < 0) && wtd_change_W > 1e-6) ){  //The change between any 2 cells can't be greater than the difference between those two cells. 
              time_step = time_step/2.;
    }
 
    else{
      arp.wtd_change_total(x,y) += ( mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W );
      params.me += mycell_change_N + mycell_change_E + mycell_change_S + mycell_change_W;
      params.N   = change_in_N_cell - arp.topo(x,y+1);
      params.S   = change_in_S_cell - arp.topo(x,y-1);
      params.W   = change_in_W_cell - arp.topo(x-1,y);
      params.E   = change_in_E_cell - arp.topo(x+1,y);

      stable = 1; 
    }


  }
return time_step;
}




void groundwater(Parameters &params, ArrayPack &arp){
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

      // Update variables with some potentially interesting values:
      //   - highest wtd
      //   - lowest wtd
      //   - greatest change in wtd during this time step
      if(arp.wtd(x,y)> max_total)
        max_total  = arp.wtd(x,y);
      else if(arp.wtd(x,y)< min_total)
        min_total  = arp.wtd(x,y);
      if(fabs(arp.wtd_change_total(x,y)) > max_change)
        max_change = fabs(arp.wtd_change_total(x,y));
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
      total_changes += arp.wtd_change_total(x,y);
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