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


const float R_v = 461.0;
const float lambda = 2.5*10E6;


int calculate_new_recharge(const Parameters &params, ArrayPack &arp){

  ofstream textfile;

  textfile.open ("text_coupled_transient.txt", std::ios_base::app);

  textfile<<"this is where we will calculate the new recharge array"<<std::endl;

  for(unsigned int i=0;i<arp.topo.size();i++){


    arp.e_sat(i) = 611 * std::exp((lambda/R_v)*((1/273.15) - (1/arp.temp(i)) ));
    arp.e_a(i) = (arp.relhum(i)/100)*arp.e_sat(i);
    if(arp.wtd(i)>0)  //if there is surface water present
      arp.evap(i) = arp.e_sat(i) - arp.e_a(i);        //TODO: get full/proper equation in here
    else              //water table is below the surface
      arp.evap(i) = arp.starting_evap(i);  //we have to reset it to not surface water conditions, else a cell that used to have surface water would still be considered that way
    //could happen in a location where climate is drying through time
  }

  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.rech(i) = arp.precip(i) - arp.evap(i);  //update recharge array with new evaporation. 
    if(arp.wtd(i)<=0)
      arp.rech(i) = std::max(arp.rech(i),0.0f);  //recharge is positive if there is no surface water


  }
  textfile<<"we have the new recharge array"<<std::endl;


//TODO: how to ensure that evaporation does not decrease the surface water to below land surface? The moment wtd dips below 0, we should turn off SW evap.
  //How to do this, particularly with long time steps where we may have e.g. 5 m water at one iteration and -5 m at the next?

textfile.close();
  return 1;
}



f2d evaporation_update(Parameters &params, ArrayPack &arp){
  richdem::Timer timer_io;



    calculate_new_recharge(params, arp);
 


  return arp.wtd;
}
