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


int TransientRun(const Parameters &params, ArrayPack &arp, const int iter, double total_changes){

  ofstream textfile;

  textfile.open (params.textfilename, std::ios_base::app);

  textfile<<"TransientRun"<<std::endl;
 
  total_changes = 0.0;                 //set values to 0
  float max_total = 0.0;
  float min_total = 0.0;
  float local_fdepth = 0.0;
  float local_kcell = 0.0;
  float local_wtd = 0.0;
  float max_change = 0.0f;
  float stability_min = 100000000000000000000.0;
  float stability_max = 0.0;

//Cycle through the entire array, adding recharge to each cell. 
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)      //skip ocean cells
      continue;

    //TODO: Does this make sense? The water table itself needs to increase by the amount of recharge, then we can get to calculating head etc. 
    //I have done this in its own for loop so that we get the new wtd everywhere first. 
    //is this also the best way to avoid recharge that is negative from evaporating too much water?

    if(arp.rech(x,y)>=0)                   //If recharge is positive, we can just add it to the wtd. 
      arp.wtd(x,y) += arp.rech(x,y);
    else{                                  //If recharge is negative, we only want to reduce wtd down to 0, and then we change the evaporative scheme to non-surface-water.
      if(arp.wtd(x,y)+arp.rech(x,y)>0)
        arp.wtd(x,y) += arp.rech(x,y);
      else{
        arp.wtd(x,y) = 0;
        arp.rech(x,y) = std::max((arp.precip(x,y) - arp.starting_evap(x,y)),0.0f); //wtd has dipped down to the surface, no more surface water in this cell. Thoughts on updating this vs leaving it until the next evaporation update time?
      }
    }
  }

//cycle through the entire array, calculating how much the water table changes in each cell per iteration. 
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)          //skip ocean cells
      continue;

    const auto my_head = arp.topo(x,y) + arp.wtd(x,y);        //just elevation head - topography plus the water table depth (negative if water table is below earth surface)               
    const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1);              
    const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1);              
    const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y);              
    const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y);              

    const auto my_kcell = kcell(x,  y,   arp);                //Get the hydraulic conductivity for our cells of interest
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);


    arp.stability_time_seconds(x,y) = 0.5 * (arp.cell_area[y]*arp.cell_area[y]) / my_kcell;  //Von Neumann stability analysis to get stable timestep
    //TODO: implement stable timestep, or remove this. 

    if(arp.stability_time_seconds(x,y)<stability_min){
      stability_min = arp.stability_time_seconds(x,y);
      local_fdepth = arp.fdepth(x,y);
      local_kcell = my_kcell;
      local_wtd = arp.wtd(x,y);
    }
    if(arp.stability_time_seconds(x,y)>stability_max)
      stability_max = arp.stability_time_seconds(x,y);


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
    QN = ((kcellN+my_kcell)/2.)*(headN-my_head)*params.deltat/params.cellsize_n_s_metres;
    QS = ((kcellS+my_kcell)/2.)*(headS-my_head)*params.deltat/params.cellsize_n_s_metres ;
    QE = ((kcellW+my_kcell)/2.)*(headW-my_head)*params.deltat/arp.cellsize_e_w_metres[y];
    QW = ((kcellE+my_kcell)/2.)*(headE-my_head)*params.deltat/arp.cellsize_e_w_metres[y] ;

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
    arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);   //update the whole wtd array at once. 
    total_changes += arp.wtd_change_total(x,y);
  }


  textfile<<"total GW changes were "<<total_changes<<std::endl;
  textfile<<"max wtd was "<<max_total<<" and min wtd was "<<min_total<<std::endl;
//  textfile<<"stability_min "<<stability_min<<" stability_max "<<stability_max<<std::endl;
  textfile<<"max GW change was "<< max_change<< std::endl;
  textfile.close();


  return 1;
}



void groundwater(Parameters &params, ArrayPack &arp){
  richdem::Timer timer_io;

  ///////////////////////////////
  //Execution Section

  int iter                  = 0;                           //Number of iterations made
  double total_changes      = 0.0;
  //Start of the main loop with condition - either number of iterations or equilibrium condition. TODO: determine a good equilibrium condition
  while(true){

    if(iter>=params.maxiter)// || abs(total_changes) < 20000.0)
      break;

    std::ofstream textfile;
    textfile.open (params.textfilename, std::ios_base::app);
  //  textfile<<"Iteration #: "<<iter<<std::endl;
    textfile.close();

    TransientRun(params, arp, iter,total_changes);
 
    iter++;
  }


}
