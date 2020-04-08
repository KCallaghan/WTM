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
           integreation of the hydraulic conductivity over flow depth.
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
  float max_total  = 0.;
  float min_total  = 0.;
  float max_change = 0.;
  
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

      //Total change in wtd for our target cell in this iteration
      arp.wtd_change_total(x,y) = ( wtd_change_N + wtd_change_S \
                                    + wtd_change_E + wtd_change_W );

      // Update variables with some potentially interesting values:
      //   - highest wtd
      //   - lowest wtd
      //   - greatest change in wtd during this time step
      if(arp.wtd(x,y)> max_total)
        max_total = arp.wtd(x,y);
      else if(arp.wtd(x,y)< min_total)
        min_total =arp.wtd(x,y);
      if(fabs(arp.wtd_change_total(x,y)) > max_change)
        max_change =fabs(arp.wtd_change_total(x,y));
    }
  }


  ////////////////
  // UPDATE WTD // 
  ////////////////

  for(int y=1;y<params.ncells_y-1;y++){
    for(int x=1;x<params.ncells_x-1; x++){
      // Update the whole wtd array at once. 
      // This is the new water table after groundwater has moved 
      // for delta_t seconds. 
      arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);   
      total_changes += arp.wtd_change_total(x,y);
    }
  }

  // Write status to text file
  textfile << "total GW changes were " << total_changes << std::endl;
  textfile << "max wtd was " << max_total << " and min wtd was " \
           << min_total << std::endl;
  textfile << "max GW change was " << max_change << std::endl;
  textfile.close();
}
