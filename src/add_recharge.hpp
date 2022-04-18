#pragma once

#include "ArrayPack.hpp"
#include "parameters.hpp"

// Method for adding recharge where porosity does not change with depth:
inline double add_recharge(const double my_rech, double my_wtd, const double my_porosity) {
  if (my_wtd >= 0) {  // all the recharge will occur above the land surface; don't worry about porosity.
    my_wtd += my_rech;
    if (my_wtd < 0) {
      my_wtd = 0;  // however, if we had evaporation, don't evaporate more surface water than is available.
    }
  } else if (my_rech > 0) {  // at least some of the water will be added into the ground, so we need to think about
                             // porosity. If my_rech was < 0 we don't add it here since there is no surface
                             // water available to evaporate.
    const double GW_space = -my_wtd * my_porosity;  // if wtd is negative, this is the amount of above-ground equivalent
                                                    // recharge that can be accommodated below ground.
    if (GW_space > my_rech) {  // all of the recharge will be below the ground; GW will not be completely filled.
      my_wtd += my_rech / my_porosity;
    } else {                        // we will have some below-ground and some above-ground water added.
      my_wtd = my_rech - GW_space;  // some is used up in GW space and the remainder goes above the ground.
    }
  }
  return my_wtd;
}
