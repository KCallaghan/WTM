#ifndef _array_pack_
#define _array_pack_

#include <richdem/common/Array2D.hpp>

namespace rd = richdem;

typedef richdem::Array2D<float>  f2d;
typedef richdem::Array2D<double>  d2d;

typedef richdem::Array2D<uint8_t>  ui82d;
typedef std::vector<double> dvec;
typedef int32_t dh_label_t;
typedef uint32_t flat_c_idx;

class ArrayPack {
 public:

  f2d porosity;
  f2d effective_storativity;
  f2d land_mask;
  f2d ice_mask;
  d2d nope;
  d2d initial_T;
  d2d scalar_array_x;
  d2d scalar_array_y;
  d2d scalar_array_x_half;
  d2d scalar_array_y_half;

  f2d ksat;
  f2d vert_ksat;

  d2d head;
  d2d transmissivity;
  d2d temp_T;

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
  d2d wtd_changed;
  d2d wtd_T;
  d2d wtd_T_iteration;
  d2d my_last_wtd;
  d2d my_prev_wtd;
  d2d original_wtd;

  dh_label_t flowdir_t;

  dvec latitude_radians;
  dvec cell_area;
  dvec cellsize_e_w_metres;
  dvec cellsize_e_w_metres_N;
  dvec cellsize_e_w_metres_S;

  rd::Array2D<dh_label_t> label; //No cells are part of a depression
  rd::Array2D<dh_label_t> final_label; //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs; //No cells flow anywhere


  void check() const;
};

#endif
