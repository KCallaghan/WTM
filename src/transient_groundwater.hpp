#ifndef transient_groundwater_h
#define transient_groundwater_h

#include "ArrayPack.hpp"
#include "parameters.hpp"

#include <richdem/common/Array2D.hpp>

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;

namespace FanDarcyGroundwater {

void update(Parameters &params, ArrayPack &arp);

}

#endif
