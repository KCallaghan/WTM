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


//Mini-function that gives the kcell, which changes through time as the water table depth changes:
///@param x         The x-coordinate of the cell in question
///
///@param y         The y-coordinate of the cell in question
///
///@param ArrayPack Global arrays. Here we use: 
///        - fdepth: The e-folding depth based on slope and temperature. 
///                  This describes the decay of kcell with depth. 
///        - wtd:    The water table depth. We use a different calculation
///                  For water tables above vs below 1.5 m below land surface.
///        - ksat:   Hydraulic conductivity, based on soil types
///@return  The kcell value for the cell in question. This is the 
///         integreation of the hydraulic conductivity over flow depth.
double kcell(const int x, const int y, const ArrayPack &arp){
  if(arp.fdepth(x,y)>0){
    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      return arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper
    else if(arp.wtd(x,y) > 0)
      return arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y));                                        //If wtd is greater than 0, max out rate of groundwater movement as though wtd were 0. The surface water will get to move in FillSpillMerge.
    else
      return arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));                             //Equation S4 from the Fan paper
  }else
    return 0;
}

///@param params   Global paramaters - we use the texfilename, run type, 
///                number of cells in the x and y directions (ncells_x and ncells_y),
///                delta_t (number of seconds in a time step), 
///                and cellsize_n_s_metres (size of a cell in the north-south direction)
///
///@param arp      Global arrays - we access land_mask, topo, wtd, wtd_change_total,
///                cellsize_e_w_metres, and cell_area.
///                land_mask is a binary representation of where land is vs where ocean is.
///                topo is the input topography, i.e. land elevation above sea level in each cell.
///                wtd is the water table depth.  
///                wtd_change_total is the amount by which wtd will change in this iteration
///                as a result of groundwater movement.
///                cellsize_e_w_metres is the size of a cell in the east-west direction.
///                cell_area is the area of the cell, needed because different cells have 
///                different areas and therefore can accommodate different water volumes.
///
///@return  An updated wtd that represents the status of the water table after groundwater
///         has been able to flow for the amount of time represented by delta_t. 
void groundwater(const Parameters &params, ArrayPack &arp){

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);
  
  textfile<<"Groundwater"<<std::endl;

  double total_changes    = 0.0;                 //reset values to 0
  float max_total  = 0.0;
  float min_total  = 0.0;
  float max_change = 0.0f;

//cycle through the entire array, calculating how much the water table changes in each cell per iteration. 
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)          //skip ocean cells
      continue;

    const auto my_head = arp.topo(x,y)   + arp.wtd(x,y);        //just elevation head - topography plus the water table depth (negative if water table is below earth surface)               
    const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1);      //and the heads for each of my neighbour cells
    const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1);              
    const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y);              
    const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y);              

    const auto my_kcell = kcell(x,  y,   arp);                //Get the hydraulic conductivity for our cells of interest
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);

    double QN = 0;             //reset all of these values to 0 before doing the calculation
    double QS = 0;
    double QE = 0;
    double QW = 0;

    double wtd_change_N = 0;
    double wtd_change_S = 0;
    double wtd_change_E = 0;
    double wtd_change_W = 0;

    arp.wtd_change_total(x,y) = 0;


//TODO: consider creating staggered grid with head gradients & kcell averages, and then doing cell-by-cell q discharges. May be faster as we don't calculate each gradient twice. 


//discharge per unit area in each direction:
    //Average hydraulic conductivity of the two cells * head difference between the two * time step (number of seconds we are moving water for, since kcell is in m/s) / cellsize (distance water is travelling)
    QN = ((kcellN+my_kcell)/2.)*(headN-my_head)*params.deltat*arp.cellsize_e_w_metres[y]/params.cellsize_n_s_metres;
    QS = ((kcellS+my_kcell)/2.)*(headS-my_head)*params.deltat*arp.cellsize_e_w_metres[y]/params.cellsize_n_s_metres ;
    QE = ((kcellW+my_kcell)/2.)*(headW-my_head)*params.deltat*params.cellsize_n_s_metres/arp.cellsize_e_w_metres[y];
    QW = ((kcellE+my_kcell)/2.)*(headE-my_head)*params.deltat*params.cellsize_n_s_metres/arp.cellsize_e_w_metres[y] ;

//divide by area of cell into which it is flowing: 
    wtd_change_N  = QN / arp.cell_area[y+1];
    wtd_change_S  = QS / arp.cell_area[y-1];
    wtd_change_E  = QE / arp.cell_area[y];
    wtd_change_W  = QW / arp.cell_area[y];

 
    arp.wtd_change_total(x,y) = (wtd_change_N + wtd_change_S + wtd_change_E + wtd_change_W);  //Total change in wtd for our target cell in this iteration


    if(arp.wtd(x,y)> max_total)         //Just some potentially interesting values - highest and lowest wtd, and greatest change in wtd this iteration. 
      max_total = arp.wtd(x,y);
    else if(arp.wtd(x,y)< min_total)
      min_total =arp.wtd(x,y);

    if(fabs(arp.wtd_change_total(x,y)) > max_change)
      max_change =fabs(arp.wtd_change_total(x,y));
  }


 for(int y=1;y<params.ncells_y-1;y++)
 for(int x=1;x<params.ncells_x-1; x++){
    arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);   //update the whole wtd array at once. This is the new water table after groundwater has moved for delta_t seconds. 
    total_changes += arp.wtd_change_total(x,y);
  }

  textfile<<"total GW changes were "<<total_changes<<std::endl;
  textfile<<"max wtd was "<<max_total<<" and min wtd was "<<min_total<<std::endl;
  textfile<<"max GW change was "<< max_change<< std::endl;
  textfile.close();
}