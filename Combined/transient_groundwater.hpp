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



typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;


//Mini-function that gives the kcell, which changes through time as the water table depth changes:
double kcell(const int x, const int y, const ArrayPack &arp){
  if(arp.fdepth(x,y)>0){
    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      return arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper
    else if(arp.wtd(x,y) > 0)
      return arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y)); 
    else
      return arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));                             //Equation S4 from the Fan paper
  } else {
    return 0;
  }
}


int TransientRun(const Parameters &params, ArrayPack &arp, const int iter, double total_changes){

  std::cout<<"TransientRun"<<std::endl;
 
  total_changes = 0.0;
  float max_total = 0.0;
  float min_total = 0.0;
  float local_fdepth = 0.0;
  float local_kcell = 0.0;
  float local_wtd = 0.0;
  float max_change = 0.0f;

  float stability_min = 100000000000000000000.0;
  float stability_max = 0.0;


//cycle through the entire array, calculating how much the water table changes in each cell per iteration. 
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)
      continue;


    const auto my_head = arp.topo(x,y) + arp.wtd(x,y) + arp.rech(x,y);                    
    const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1) + arp.rech(x,y+1);              
    const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1) + arp.rech(x,y-1);              
    const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y) + arp.rech(x-1,y);              
    const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y) + arp.rech(x+1,y);              

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);


    arp.stability_time_seconds(x,y) = 0.5 * (arp.cell_area[y]*arp.cell_area[y]) / my_kcell;  //Von Neumann stability analysis to get stable timestep
    //TODO: implement stable timestep


    if(arp.stability_time_seconds(x,y)<stability_min){
      stability_min = arp.stability_time_seconds(x,y);
      local_fdepth = arp.fdepth(x,y);
      local_kcell = my_kcell;
      local_wtd = arp.wtd(x,y);
    }
    if(arp.stability_time_seconds(x,y)>stability_max)
      stability_max = arp.stability_time_seconds(x,y);



    double QN = 0;
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
    QN = ((kcellN+my_kcell)/2.)*(headN-my_head)*params.deltat/params.cellsize_n_s_metres;
    QS = ((kcellS+my_kcell)/2.)*(headS-my_head)*params.deltat/params.cellsize_n_s_metres ;
    QE = ((kcellW+my_kcell)/2.)*(headW-my_head)*params.deltat/arp.cellsize_e_w_metres[y];
    QW = ((kcellE+my_kcell)/2.)*(headE-my_head)*params.deltat/arp.cellsize_e_w_metres[y] ;

//multiply by number of seconds we are moving water for, divide by distance water is travelling, divide by area of cell into which it is flowing: 
    wtd_change_N  = QN / arp.cell_area[y+1];
    wtd_change_S  = QS / arp.cell_area[y-1];
    wtd_change_E  = QE / arp.cell_area[y];
    wtd_change_W  = QW / arp.cell_area[y];

 
    arp.wtd_change_total(x,y) = (wtd_change_N + wtd_change_S + wtd_change_E + wtd_change_W);


    if(arp.wtd(x,y)> max_total)
      max_total = arp.wtd(x,y);
    else if(arp.wtd(x,y)< min_total)
      min_total =arp.wtd(x,y);

    if(fabs(arp.wtd_change_total(x,y)) > max_change){
      max_change =fabs(arp.wtd_change_total(x,y));
    }
  


  }


 for(int y=1;y<params.ncells_y-1;y++)
 for(int x=1;x<params.ncells_x-1; x++){
    if(arp.ksat(x,y) == 0)
      continue;

    arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);   
 
    total_changes += arp.wtd_change_total(x,y);

  }

std::cout<<"total changes were "<<total_changes<<std::endl;
std::cout<<"max wtd was "<<max_total<<" and min wtd was "<<min_total<<std::endl;
std::cout<<"stability_min "<<stability_min<<" stability_max "<<stability_max<<std::endl;
//std::cout<<"min kcell "<<local_kcell<<" min wtd "<<local_wtd<<" min fdepth "<<local_fdepth<<std::endl;

std::cout<<"max change was "<< max_change<< std::endl;


  return 1;
}



f2d transient(Parameters &params, ArrayPack &arp){
  richdem::Timer timer_io;



  ///////////////////////////////
  //Execution Section

  int iter                  = 0;                           //Number of iterations made
  double total_changes      = 0.0;
  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium.
  while(true){

    if(iter>=params.maxiter)// || abs(total_changes) < 20000.0)
      break;

    std::cerr<<"Iteration #: "<<iter<<std::endl;

    TransientRun(params, arp, iter,total_changes);
 
    iter++;
  }


  return arp.wtd;
}
