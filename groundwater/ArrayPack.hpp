#ifndef _array_pack_
#define _array_pack_

#include <richdem/common/Array2D.hpp>
#include <vector>

typedef richdem::Array2D<float>  f2d;
typedef std::vector<double> dvec;

class ArrayPack {
 public:
  f2d ksat;   
  f2d wtd;
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

  richdem::Array2D<bool> done_new;   //Indicates which cells must still be processed
  richdem::Array2D<bool> done_old;   //Indicates which cells must still be processed
  richdem::Array2D<bool> land;

  dvec xlat;
  dvec alpha;
  dvec alphamonth;

  void check() const;
};

#endif
