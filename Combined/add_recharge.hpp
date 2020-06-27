#include "ArrayPack.hpp"
#include "parameters.hpp"


void add_recharge(const int x, const int y, Parameters &params, ArrayPack &arp){

  double volume_change = arp.rech(x,y)/31536000. * params.deltat * arp.cell_area[y];
  double GW_portion = 0;


  if(arp.wtd(x,y) >= 0){  //all the recharge will occur above the land surface; don't worry about porosity. 
    arp.wtd(x,y) += arp.rech(x,y)/31536000. * params.deltat;
    if(arp.wtd(x,y) < 0)
      arp.wtd(x,y) = 0; //however, if we had evaporation, don't evaporate more surface water than is available. 
  }
  else if(volume_change>0){ //at least some of the water will be added into the ground, so we need to think about porosity. 
    arp.wtd(x,y) = arp.fdepth(x,y) * log(exp(arp.wtd(x,y) / arp.fdepth(x,y)) \
          + volume_change / ( arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)) );

    if(arp.wtd(x,y) > 0 || (exp(arp.wtd(x,y)/arp.fdepth(x,y)) <  (volume_change / (arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)))  ) ){ //part of the change is below the ground, part above. So we need to calculate how much of the water was used up in the ground, i.e. the portion between the starting wtd and 0.
      GW_portion = -arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y) * \
                       (exp(arp.wtd(x,y) / arp.fdepth(x,y)) - 1);
                       //this is the volume of water used up in filling in the ground. 
                       //volume_change - GW_portion is left over to fill surface water. 
      arp.wtd(x,y) = ((volume_change - GW_portion) / arp.cell_area[y]);


    }
  }


}



  