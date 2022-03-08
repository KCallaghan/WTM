#include "ArrayPack.hpp"
#include "parameters.hpp"

//Method for adding recharge where porosity does not change with depth:
double add_recharge(const double deltat, const double my_rech, double my_wtd, const int my_mask, const double my_porosity){
  constexpr double seconds_in_a_year = 31536000.;

  double rech_change = my_rech/seconds_in_a_year * deltat;

  if(my_wtd >= 0 && my_mask == 1){  //all the recharge will occur above the land surface; don't worry about porosity.
    my_wtd += rech_change;
    if(my_wtd < 0)
      my_wtd = 0; //however, if we had evaporation, don't evaporate more surface water than is available.
  }
  else if(rech_change>0 && my_mask == 1){ //at least some of the water will be added into the ground, so we need to think about porosity. If rech_change was < 0 we don't add it here since there is no surface water available to evaporate.
    double GW_space = -my_wtd * my_porosity;  //if wtd is negative, this is the amount of above-ground equivalent recharge that can be accommodated below ground.
    if(GW_space > rech_change) //all of the recharge will be below the ground; GW will not be completely filled.
      my_wtd += rech_change / my_porosity;
    else  //we will have some below-ground and some above-ground water added.
      my_wtd = rech_change - GW_space; //some is used up in GW space and the remainder goes above the ground.
  }
  return my_wtd;
}
