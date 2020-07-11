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
///We use Dalton's Law to calculate evaporation
///in cells that contain surface water. 
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

  #pragma omp parallel for 
  for(unsigned int i=0;i<arp.topo.size();i++){

    arp.e_sat(i) = 0.611 * std::exp((17.3*arp.ground_temp(i))\
    /(arp.ground_temp(i)+237.3));

    arp.e_a(i) = arp.relhum(i) * 0.611 * std::exp((17.3*arp.temp(i))\
    /(arp.temp(i)+237.3));

    arp.surface_evap(i) = (K_e*arp.wind_speed(i))*(arp.e_sat(i) - arp.e_a(i)) * 31536000;        //convert from m/s to m/yr        

    if(arp.wtd(i)>0)  //if there is surface water present
      arp.rech(i) = arp.precip(i) - arp.surface_evap(i);
    else{              //water table is below the surface
      arp.rech(i) = arp.precip(i) - arp.starting_evap(i);
      if(arp.rech(i) <0)    //Recharge is always positive. 
        arp.rech(i) = 0.0f;
     }
      //we have to reset it to not surface water conditions, 
    //else a cell that used to have surface water would still be 
    //considered that way.
   //could happen in a location where climate is drying through time
  }
}