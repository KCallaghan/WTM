#include "ArrayPack.hpp"
#include "parameters.hpp"

//Method for adding recharge where porosity does not change with depth:
double add_recharge(const double deltat, const double my_rech, double my_wtd, const int my_mask, const double my_porosity,Parameters &params,const double cell_area){

  double rech_change = my_rech/31536000. * deltat;

  if(my_wtd >= 0 && my_mask == 1){  //all the recharge will occur above the land surface; don't worry about porosity.
    my_wtd += rech_change;
    params.total_added_recharge += rech_change*cell_area;
    if(my_wtd < 0){
      params.total_added_recharge += my_wtd*cell_area;
      my_wtd = 0; //however, if we had evaporation, don't evaporate more surface water than is available.
    }
  }
  else if(rech_change>0 && my_mask == 1){ //at least some of the water will be added into the ground, so we need to think about porosity. If rech_change was < 0 we don't add it here since there is no surface water available to evaporate.
    params.total_added_recharge += rech_change*cell_area;
    double GW_space = -my_wtd * my_porosity;  //if wtd is negative, this is the amount of above-ground equivalent recharge that can be accommodated below ground.
    if(GW_space > rech_change) //all of the recharge will be below the ground; GW will not be completely filled.
      my_wtd += rech_change / my_porosity;
    else  //we will have some below-ground and some above-ground water added.
      my_wtd = rech_change - GW_space; //some is used up in GW space and the remainder goes above the ground.
  }
  return my_wtd;
}




//Method for adding recharge if depth-variable porosity is used:
//void add_recharge(const int x, const int y, double deltat, ArrayPack &arp){
//
//  double volume_change = arp.rech(x,y)/31536000. * deltat * arp.cell_area[y];
//  double GW_portion = 0;
//
//
//  if(arp.wtd(x,y) >= 0){  //all the recharge will occur above the land surface; don't worry about porosity.
//    arp.wtd(x,y) += arp.rech(x,y)/31536000. * deltat;
//    if(arp.wtd(x,y) < 0)
//      arp.wtd(x,y) = 0; //however, if we had evaporation, don't evaporate more surface water than is available.
//  }
//  else if(volume_change>0){ //at least some of the water will be added into the ground, so we need to think about porosity.
//    arp.wtd(x,y) = arp.fdepth(x,y) * log(exp(arp.wtd(x,y) / arp.fdepth(x,y))
//          + volume_change / ( arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)) );
//
//    if(arp.wtd(x,y) > 0 || (exp(arp.wtd(x,y)/arp.fdepth(x,y)) <  (volume_change / (arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y)))  ) ){ //part of the change is below the ground, part above. So we need to calculate how much of the water was used up in the ground, i.e. the portion between the starting wtd and 0.
//      GW_portion = -arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y) *
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