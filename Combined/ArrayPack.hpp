#ifndef _array_pack_
#define _array_pack_

#include <richdem/common/Array2D.hpp>

namespace rd = richdem;

typedef richdem::Array2D<float>  f2d;
typedef std::vector<double> dvec;
typedef int32_t dh_label_t;

class ArrayPack {
 public:
  f2d ksat;  
  f2d vert_ksat;  
  f2d porosity;  
  f2d slope_start;
  f2d slope_end;
  f2d fdepth_start; 
  f2d fdepth_end;
  f2d precip_start; 
  f2d precip_end;
  f2d temp_start; 
  f2d temp_end;
  f2d topo_start;   
  f2d winter_temp_start; 
  f2d winter_temp_end;
  f2d winter_temp;

  f2d topo_end;
  f2d ground_temp_start;
  f2d ground_temp_end;
  f2d ground_temp;
  f2d wind_speed_start;
  f2d wind_speed_end;
  f2d wind_speed;
  f2d starting_evap_start; 
  f2d starting_evap_end;
  f2d relhum_start; 
  f2d relhum_end;
  f2d wtd;
  f2d infiltration_array;
  f2d surface_array;
  f2d starting_rech;

  f2d fdepth;
  f2d precip;
  f2d temp;
  f2d topo;
  f2d starting_evap;
  f2d relhum;
  f2d slope;
  f2d land_mask;

  f2d wtd_old;
  f2d wtd_mid;
  f2d rech;
  f2d runoff;
  f2d head;
  f2d kcell;
  f2d evap;
  f2d e_sat;
  f2d e_a;
  f2d surface_evap;
  f2d wtd_change_total;

  f2d stability_time_seconds;

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
