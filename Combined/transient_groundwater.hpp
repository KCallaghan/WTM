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

void groundwater(const Parameters &params, ArrayPack &arp){
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

  // Declare variables updated in the loop
  double wtd_change_N;
  double wtd_change_S;
  double wtd_change_E;
  double wtd_change_W;

  // Declare status variables and set initial values to 0
  double total_changes = 0.;
  float max_total      = 0.;
  float min_total      = 0.;
  float max_change     = 0.;

  float stability_min = 100000000000000000000.0;
  float stability_max = 0.0;

  
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

      // Elevation head - topography plus the water table depth (negative if 
      // water table is below earth surface)               
      const auto my_head = arp.topo(x,y) + arp.wtd(x,y);
      // heads for each of my neighbour cells
      const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1);
      const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1);              
      const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y);              
      const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y);              

      // Get the hydraulic conductivity for our cells of interest
      const auto kN = ( kcell(x,  y,   arp) + kcell(x,  y+1, arp) ) / 2.;
      const auto kS = ( kcell(x,  y,   arp) + kcell(x,  y-1, arp) ) / 2.;                
      const auto kW = ( kcell(x,  y,   arp) + kcell(x-1,y,   arp) ) / 2.;
      const auto kE = ( kcell(x,  y,   arp) + kcell(x+1,y,   arp) ) / 2.;

      float water_height_N = 0.0;
      float water_height_S = 0.0;
      float water_height_W = 0.0;
      float water_height_E = 0.0;




          arp.stability_time_seconds(x,y) = 0.5 * (arp.cell_area[y]*arp.cell_area[y]) / kcell(x,y,arp);  //Von Neumann stability analysis to get stable timestep
    //TODO: implement stable timestep, or remove this. 

      if(arp.stability_time_seconds(x,y)<stability_min){
        stability_min = arp.stability_time_seconds(x,y);
      }
      if(arp.stability_time_seconds(x,y)>stability_max)
        stability_max = arp.stability_time_seconds(x,y);




      // Change in water-table depth.
      // (1) Discharge across cell boundaries
      // Average hydraulic conductivity of the two cells * 
      // head difference between the two / distance (i.e., dH/dx_i) *
      // width of cell across which the water is discharged *
      // time step
      // (2) Divide by the area of the given cell: maps water volume 
      // increases/decreases to change in head
      wtd_change_N = kN * (headN - my_head) / params.cellsize_n_s_metres \
                      * arp.cellsize_e_w_metres[y] * params.deltat \
                      / arp.cell_area[y];
      wtd_change_S = kS * (headS - my_head) / params.cellsize_n_s_metres \
                      * arp.cellsize_e_w_metres[y] * params.deltat \
                      / arp.cell_area[y];
      wtd_change_E = kE * (headE - my_head) / arp.cellsize_e_w_metres[y] \
                      * params.cellsize_n_s_metres * params.deltat \
                      / arp.cell_area[y];
      wtd_change_W = kW * (headW - my_head) / arp.cellsize_e_w_metres[y] \
                      * params.cellsize_n_s_metres * params.deltat \
                      / arp.cell_area[y];

      //Now think about porosity for each of these directions:
      if(wtd_change_N > 0){  //the current cell is receiving water from the North
        if(arp.wtd(x,y+1) > 0){
          water_height_N = wtd_change_N; //if it is only surface water, it does not change. 
          if(wtd_change_N > arp.wtd(x,y+1)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_N = arp.wtd(x,y+1) + (wtd_change_N - arp.wtd(x,y+1))*arp.porosity(x,y+1); 
          }    
        }
        else if(arp.wtd(x,y+1) < 0){ //all groundwater that needs porosity
          water_height_N = wtd_change_N*arp.porosity(x,y+1);
        }
      }
      else{ //the current cell is giving water to the North. 
        if(arp.wtd(x,y) > 0){
          water_height_N = wtd_change_N; //if it is only surface water, it does not change. 
          if(fabs(wtd_change_N) > arp.wtd(x,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_N = -arp.wtd(x,y) - (-wtd_change_N - arp.wtd(x,y))*arp.porosity(x,y); 
          }    
        }
        else if(arp.wtd(x,y) < 0){ //all groundwater that needs porosity
          water_height_N = wtd_change_N*arp.porosity(x,y);
        }
      }





      if(wtd_change_S > 0){  //the current cell is receiving water from the South
        if(arp.wtd(x,y-1) > 0){
          water_height_S = wtd_change_S; //if it is only surface water, it does not change. 
          if(wtd_change_S > arp.wtd(x,y-1)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_S = arp.wtd(x,y-1) + (wtd_change_S - arp.wtd(x,y-1))*arp.porosity(x,y-1); 
          }    
        }
        else if(arp.wtd(x,y-1) < 0){ //all groundwater that needs porosity
          water_height_S = wtd_change_S*arp.porosity(x,y-1);
        }
      }
      else{ //the current cell is giving water to the South. 
        if(arp.wtd(x,y) > 0){
          water_height_S = wtd_change_S; //if it is only surface water, it does not change. 
          if(fabs(wtd_change_S) > arp.wtd(x,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_S = -arp.wtd(x,y) - (-wtd_change_S - arp.wtd(x,y))*arp.porosity(x,y); 
          }    
        }
        else if(arp.wtd(x,y) < 0){ //all groundwater that needs porosity
          water_height_S = wtd_change_S*arp.porosity(x,y);
        }
      }



      if(wtd_change_W > 0){  //the current cell is receiving water from the West
        if(arp.wtd(x-1,y) > 0){
          water_height_W = wtd_change_W; //if it is only surface water, it does not change. 
          if(wtd_change_W > arp.wtd(x-1,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_W = arp.wtd(x-1,y) + (wtd_change_W - arp.wtd(x-1,y))*arp.porosity(x-1,y); 
          }    
        }
        else if(arp.wtd(x-1,y) < 0){ //all groundwater that needs porosity
          water_height_W = wtd_change_W*arp.porosity(x-1,y);
        }
      }
      else{ //the current cell is giving water to the West. 
        if(arp.wtd(x,y) > 0){
          water_height_W = wtd_change_W; //if it is only surface water, it does not change. 
          if(fabs(wtd_change_W) > arp.wtd(x,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_W = -arp.wtd(x,y) - (-wtd_change_W - arp.wtd(x,y))*arp.porosity(x,y); 
          }    
        }
        else if(arp.wtd(x,y) < 0){ //all groundwater that needs porosity
          water_height_W = wtd_change_W*arp.porosity(x,y);
        }
      }



      if(wtd_change_E > 0){  //the current cell is receiving water from the East
        if(arp.wtd(x+1,y) > 0){
          water_height_E = wtd_change_E; //if it is only surface water, it does not change. 
          if(wtd_change_E > arp.wtd(x+1,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_E = arp.wtd(x+1,y) + (wtd_change_E - arp.wtd(x+1,y))*arp.porosity(x+1,y); 
          }    
        }
        else if(arp.wtd(x+1,y) < 0){ //all groundwater that needs porosity
          water_height_E = wtd_change_E*arp.porosity(x+1,y);
        }
      }
      else{ //the current cell is giving water to the East. 
        if(arp.wtd(x,y) > 0){
          water_height_E = wtd_change_E; //if it is only surface water, it does not change. 
          if(fabs(wtd_change_E) > arp.wtd(x,y)){ 
          //some surface water that does not think about porosity, and some groundwater that does. 
            water_height_E = -arp.wtd(x,y) - (-wtd_change_E - arp.wtd(x,y))*arp.porosity(x,y); 
          }    
        }
        else if(arp.wtd(x,y) < 0){ //all groundwater that needs porosity
          water_height_E = wtd_change_E*arp.porosity(x,y);
        }
      }



  //    arp.wtd_change_total(x,y) = ( wtd_change_N + wtd_change_E \
                                    + wtd_change_S + wtd_change_W );
//std::cout<<"wtd_changes: N "<<wtd_change_N<<" S "<<wtd_change_S<<" E "<<wtd_change_E<<" W "<<wtd_change_W<<std::endl;
//std::cout<<"water_heights: N "<<water_height_N<<" S "<<water_height_S<<" E "<<water_height_E<<" W "<<water_height_W<<std::endl;

      //Total change in wtd for our target cell in this iteration
      arp.wtd_change_total(x,y) = ( water_height_N + water_height_E \
                                    + water_height_S + water_height_W );

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

  std::cout<<"stability_min "<<stability_min<<" stability_max "<<stability_max<<std::endl;

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


      //now we must once again think about porosity. 
      if(arp.wtd(x,y)>0){
        if(arp.wtd(x,y) + arp.wtd_change_total(x,y) > 0){
          //it's all surface water, so no multiplier needed
          arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);  
        }
        else{ //some surface and some groundwater, and we are losing water. 
          arp.wtd_change_total(x,y) += arp.wtd(x,y); //remove the surface water part. wtd_change_total is negative, so add the positive wtd to see how much is left.
          arp.wtd(x,y) = arp.wtd_change_total(x,y) / arp.porosity(x,y); // wtd_change_total was negative, now so is wtd. 

        } 
      }
      else{ //the starting wtd is negative. 
        if(arp.wtd(x,y) + arp.wtd_change_total(x,y) < 0){
          //it is staying negative, so all is affected by porosity. 
          arp.wtd(x,y) += arp.wtd_change_total(x,y) / arp.porosity(x,y);
        }
       else{ //water table is getting filled in and some surface water added. 
          arp.wtd_change_total(x,y) += arp.wtd(x,y)*arp.porosity(x,y); //remove the surface water part. wtd_change_total is positive, so add the negative wtd to see how much is left.
          arp.wtd(x,y) = arp.wtd_change_total(x,y); // wtd_change_total was positive, now so is wtd. 

        }

      }
  //  std::cout<<"wtd "<<arp.wtd(x,y)<<" change "<<arp.wtd_change_total(x,y)<<" porosity "<<arp.porosity(x,y)<<" x "<<x<<" y "<<y<<std::endl;
    }
  }

  // Write status to text file
  textfile << "total GW changes were " << total_changes << std::endl;
  textfile << "max wtd was " << max_total << " and min wtd was " \
           << min_total << std::endl;
  textfile << "max GW change was " << max_change << std::endl;
  textfile.close();
}
