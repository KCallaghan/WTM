#pragma once

#include <richdem/common/Array2D.hpp>
#include "dephier.hpp"

typedef richdem::Array2D<float> f2d;
typedef richdem::Array2D<double> d2d;

typedef std::vector<double> dvec;

struct ArrayPack {
  // input data files - transient:

  f2d evap_start;
  f2d evap_end;
  f2d open_water_evap_start;
  f2d open_water_evap_end;
  f2d precip_start;
  f2d precip_end;
  f2d runoff_ratio_start;
  f2d runoff_ratio_end;
  f2d slope_start;
  f2d slope_end;
  f2d topo_start;
  f2d topo_end;
  f2d winter_temp_start;
  f2d winter_temp_end;

  // input data files - equilibrium:

  f2d evap;
  f2d open_water_evap;
  f2d precip;
  f2d runoff_ratio;
  f2d slope;
  f2d topo;
  f2d winter_temp;

  // input data files - both

  f2d ksat;
  f2d land_mask;
  f2d porosity;
  f2d vert_ksat;

  // other data storage arrays:

  d2d effective_storativity;
  d2d fdepth;
  d2d fdepth_start;
  d2d fdepth_end;
  d2d infiltration_array;
  d2d rech;
  d2d runoff;
  d2d transmissivity;

  // arrays recording various states of water table depth:

  d2d wtd;
  d2d wtd_mid;
  d2d wtd_old;

  // arrays used for calculations in code:

  dvec cell_area;
  dvec cellsize_e_w_metres;

  // labels and flow directions:

  richdem::Array2D<richdem::dephier::dh_label_t> label;        // No cells are part of a depression
  richdem::Array2D<richdem::dephier::dh_label_t> final_label;  // No cells are part of a depression
  richdem::Array2D<richdem::flowdir_t> flowdirs;               // No cells flow anywhere

  // Cumulative state variables
  double total_added_recharge = 0;
  double total_loss_to_ocean  = 0;

  void check() const;
};
