#ifndef transient_groundwater_h
#define transient_groundwater_h

#include "../common/netcdf.hpp"
#include "ArrayPack.hpp"
#include "parameters.hpp"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;

namespace FanDarcyGroundwater {

void update(const Parameters &params, ArrayPack &arp);

}

#endif
