#ifndef _array_pack_
#define _array_pack_

#include <richdem/common/Array2D.hpp>
#include <vector>
#include "../common/netcdf.hpp"
#include <iostream>
#include <string>
#include <stdexcept>
//#include "dephier.hpp"

namespace rd = richdem;
//namespace dh = richdem::dephier;


typedef richdem::Array2D<float>  f2d;
typedef std::vector<double> dvec;
typedef int32_t dh_label_t;

//typedef richdem::Array2D<dh::dh_label_t> lab_arr;

class ArrayPack {
 public:
  f2d ksat;   
  f2d wtd;
  f2d wtd_old;
  f2d kcell;
  f2d head;
  f2d wtd_mid;


  f2d wtd_new;

  f2d topo_start;   
  f2d rech_start;   
  f2d fslope_start; 
  f2d temp_start;   

  f2d topo_end;
  f2d rech_end;
  f2d fslope_end;
  f2d temp_end;

  f2d topo;
  f2d rech;
  f2d fdepth;
  f2d temp;
  f2d delta;
  f2d e_sat;
  f2d e_a;
  f2d evap;
  f2d starting_evap;
  f2d relhum;
  f2d precip;
  f2d land_mask;
  f2d total;
  f2d wtd_change_total;
  f2d stability_time_seconds;


  dh_label_t label;
  dh_label_t final_label;
  dh_label_t flowdir_t;



  richdem::Array2D<bool> done_new;   //Indicates which cells must still be processed
  richdem::Array2D<bool> done_old;   //Indicates which cells must still be processed
 // richdem::Array2D<int> land;

  dvec latitude_radians;
  dvec cell_area;
  dvec cellsize_e_w_metres;
  dvec cellsize_e_w_metres_N;
  dvec cellsize_e_w_metres_S;
  dvec threshold_array;

  dvec temp_east_west;
  dvec temp_north;
  dvec temp_south;
  dvec alpha;
  dvec alphamonth;

  void check() const;
};

#endif
