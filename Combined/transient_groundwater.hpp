#ifndef transient_groundwater_h
#define transient_groundwater_h

#include "ArrayPack.hpp"
#include "parameters.hpp"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;

namespace FanDarcyGroundwater {

void update(Parameters &params, ArrayPack &arp);

}

#endif
