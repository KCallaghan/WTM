#include "ArrayPack.hpp"
#include "parameters.hpp"

//Method for adding recharge where porosity does not change with depth:
void add_recharge(const int x, const int y, Parameters &params, ArrayPack &arp){

  double rech_change = arp.rech(x,y)/31536000. * params.deltat;

  if(arp.wtd(x,y) >= 0){  //all the recharge will occur above the land surface; don't worry about porosity. 
    arp.wtd(x,y) += rech_change;
    if(arp.wtd(x,y) < 0)
      arp.wtd(x,y) = 0; //however, if we had evaporation, don't evaporate more surface water than is available. 
  }
  else if(rech_change>0){ //at least some of the water will be added into the ground, so we need to think about porosity. If rech_change was < 0 we don't add it here since there is no surface water available to evaporate.
    double GW_space = -arp.wtd(x,y) * arp.porosity(x,y);  //if wtd is negative, this is the amount of above-ground equivalent recharge that can be accommodated below ground.
    if(GW_space > rech_change){ //all of the recharge will be below the ground; GW will not be completely filled. 
      arp.wtd(x,y) += rech_change / arp.porosity(x,y);
    }
    else{  //we will have some below-ground and some above-ground water added. 
      arp.wtd(x,y) = rech_change - GW_space; //some is used up in GW space and the remainder goes above the ground. 
    }
  }
}


//Method for adding recharge if depth-variable porosity is used:
//void add_recharge(const int x, const int y, Parameters &params, ArrayPack &arp){
//
//  double volume_change = arp.rech(x,y)/31536000. * params.deltat * arp.cell_area[y];
//  double GW_portion = 0;
//
//
//  if(arp.wtd(x,y) >= 0){  //all the recharge will occur above the land surface; don't worry about porosity. 
//    arp.wtd(x,y) += arp.rech(x,y)/31536000. * params.deltat;
//    if(arp.wtd(x,y) < 0)
//      arp.wtd(x,y) = 0; //however, if we had evaporation, don't evaporate more surface water than is available. 
//  }
//  else if(volume_change>0){ //at least some of the water will be added into the ground, so we need to think about porosity. 
//    arp.wtd(x,y) = arp.fdepth(x,y) * log(exp(arp.wtd(x,y) / arp.fdepth(x,y)) \
//          + volume_change / ( arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)) );
//
//    if(arp.wtd(x,y) > 0 || (exp(arp.wtd(x,y)/arp.fdepth(x,y)) <  (volume_change / (arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)))  ) ){ //part of the change is below the ground, part above. So we need to calculate how much of the water was used up in the ground, i.e. the portion between the starting wtd and 0.
//      GW_portion = -arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y) * \
//                       (exp(arp.wtd(x,y) / arp.fdepth(x,y)) - 1);
//                       //this is the volume of water used up in filling in the ground. 
//                       //volume_change - GW_portion is left over to fill surface water. 
//      arp.wtd(x,y) = ((volume_change - GW_portion) / arp.cell_area[y]);
//    }
//  }
//}
//
//
//  