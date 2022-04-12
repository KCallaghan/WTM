#pragma once

#include <richdem/common/Array2D.hpp>
#include "dephier.hpp"

typedef richdem::Array2D<float> f2d;
typedef richdem::Array2D<double> d2d;

typedef richdem::Array2D<uint8_t> ui82d;
typedef std::vector<double> dvec;

struct ArrayPack {
  f2d porosity;
  d2d effective_storativity;
  f2d land_mask;
  f2d ice_mask;
  d2d scalar_array_x;
  d2d scalar_array_y;

  f2d ksat;
  f2d vert_ksat;

  d2d transmissivity;

  f2d topo_start;
  f2d topo_end;
  f2d topo;

  d2d runoff;
  d2d rech;
  d2d evap;
  d2d infiltration_array;

  f2d slope_start;
  f2d slope_end;
  d2d fdepth_start;
  d2d fdepth_end;
  f2d precip_start;
  f2d precip_end;
  f2d starting_evap_start;
  f2d starting_evap_end;
  f2d open_water_evap_start;
  f2d open_water_evap_end;
  f2d winter_temp_start;
  f2d winter_temp_end;

  f2d slope;
  d2d fdepth;
  f2d open_water_evap;
  f2d precip;
  f2d starting_evap;
  f2d winter_temp;

  d2d wtd;
  d2d wtd_mid;
  d2d wtd_old;
  d2d wtd_T;
  d2d wtd_T_iteration;
  d2d original_wtd;

  dvec latitude_radians;
  dvec cell_area;
  dvec cellsize_e_w_metres;

  richdem::Array2D<richdem::dephier::dh_label_t> label;        // No cells are part of a depression
  richdem::Array2D<richdem::dephier::dh_label_t> final_label;  // No cells are part of a depression
  richdem::Array2D<richdem::flowdir_t> flowdirs;               // No cells flow anywhere

  void check() const;
};
