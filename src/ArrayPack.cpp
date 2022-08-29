#include "ArrayPack.hpp"
#include <cassert>

// Check to see that all initalised arrays have the same dimensions, to help avoid incorrect input data
void ArrayPack::check() const {
  const auto compare_dimensions = [](const auto& arr1, const auto& arr2) {
    if (arr1.width() != arr2.width() || arr1.height() != arr2.height()) {
      throw std::runtime_error("Array dimensions don't match!");
    }
  };
  compare_dimensions(topo, evap);
  compare_dimensions(topo, open_water_evap);
  compare_dimensions(topo, precip);
  compare_dimensions(topo, runoff_ratio);
  compare_dimensions(topo, slope);
  compare_dimensions(topo, winter_temp);
  compare_dimensions(topo, ksat);
  compare_dimensions(topo, land_mask);
  compare_dimensions(topo, porosity);
  compare_dimensions(topo, effective_storativity);
  compare_dimensions(topo, evap);
  compare_dimensions(topo, fdepth);
  compare_dimensions(topo, infiltration_array);
  compare_dimensions(topo, rech);
  compare_dimensions(topo, runoff);
  compare_dimensions(topo, transmissivity);
  compare_dimensions(topo, wtd);
  compare_dimensions(topo, wtd_mid);
  compare_dimensions(topo, wtd_old);
  compare_dimensions(topo, label);
  compare_dimensions(topo, final_label);
  compare_dimensions(topo, flowdirs);
  if (!fdepth_end.empty()) {
    compare_dimensions(topo, evap_start);
    compare_dimensions(topo, evap_end);
    compare_dimensions(topo, open_water_evap_start);
    compare_dimensions(topo, open_water_evap_end);
    compare_dimensions(topo, precip_start);
    compare_dimensions(topo, precip_end);
    compare_dimensions(topo, runoff_ratio_start);
    compare_dimensions(topo, runoff_ratio_end);
    compare_dimensions(topo, slope_start);
    compare_dimensions(topo, slope_end);
    compare_dimensions(topo, topo_start);
    compare_dimensions(topo, topo_end);
    compare_dimensions(topo, winter_temp_start);
    compare_dimensions(topo, winter_temp_end);
    compare_dimensions(topo, fdepth_start);
    compare_dimensions(topo, fdepth_end);
  }
  if (!vert_ksat.empty()) {
    compare_dimensions(topo, vert_ksat);
  }
}
