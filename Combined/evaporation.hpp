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


///update the amount of evaporation that occurs in each cell. 
///We are calculating the vapour pressure deficit as a 
///coefficient for evaporation. 
///surface_evap is used by fill_spill_merge for evaporation 
///that occurs as water moves cell-by-cell across the land
///surface. 
///
///The amount of runoff input to a cell is also influenced
///by evaporation: in cases where surface water is present
///in a cell, the evaporation may be larger, vs cells without
///surface water. However, runoff is always a minimum of zero. 
///TODO: does this mean I am doubling up on evaporation in lakes where I should not be?
//void evaporation_update(Parameters &params, ArrayPack &arp){
//  richdem::Timer timer_io;
//
//  const float R_v = 461.0;
//  const float lambda = 2.5*10E6;
//
//  for(unsigned int i=0;i<arp.topo.size();i++){
//
//    arp.e_sat(i) = 611 * std::exp((lambda/R_v)*((1/273.15) - (1/(arp.temp(i)+273.15)) ));
//    arp.e_a(i) = arp.relhum(i)*arp.e_sat(i);
//    arp.surface_evap(i) = (arp.e_sat(i) - arp.e_a(i));        //TODO: get full/proper equation in here. Or is this equation okay?
//
//    if(arp.wtd(i)>0)  //if there is surface water present
//      arp.evap(i) = (arp.e_sat(i) - arp.e_a(i));        //TODO: get full/proper equation in here
//    else              //water table is below the surface
//      arp.evap(i) = arp.starting_evap(i);  //we have to reset it to not surface water conditions, else a cell that used to have surface water would still be considered that way
//   //could happen in a location where climate is drying through time
//  }
////
//  for(unsigned int i=0;i<arp.topo.size();i++){
//    arp.runoff(i) = (arp.precip(i)-arp.evap(i));
////  //  arp.runoff(i) = std::max(arp.runoff(i),0.0f);  //recharge is positive if there is no surface water
//  }
//}
//




void evaporation_update(Parameters &params, ArrayPack &arp){
  float atm_p = 101.3;
  float p_a = 1.220;
  float p_w = 1000;
  float k = 0.4;
  float z = 2;
  float z_d = 0;
  float z_0 = 2.3*(std::pow(10,-4));

  float log_bracket = log((z-z_d)/z_0);

  float K_e = (0.622*p_a*k*k)/(atm_p*p_w*std::pow(log_bracket,2));

  float v_a = 1;







for(unsigned int i=0;i<arp.topo.size();i++){

    arp.e_sat(i) = 0.611 * std::exp((17.3*arp.temp(i))/(arp.temp(i)+237.3));
    arp.e_a(i) = arp.relhum(i)*arp.e_sat(i);
    arp.surface_evap(i) = (K_e*v_a)*(arp.e_sat(i) - arp.e_a(i));        //TODO: get full/proper equation in here. Or is this equation okay?

    if(arp.wtd(i)>0)  //if there is surface water present
      arp.runoff(i) = arp.precip(i) - arp.surface_evap(i);
    else              //water table is below the surface
      arp.runoff(i) = arp.starting_rech(i);  //we have to reset it to not surface water conditions, else a cell that used to have surface water would still be considered that way
   //could happen in a location where climate is drying through time
  }
}