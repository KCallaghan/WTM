#include "ArrayPack.hpp"
#include <cassert>

// Check to see that all initalised arrays have the same dimensions, to help avoid incorrect input data
void ArrayPack::check() const {
  assert(topo.width() == evap.width() && topo.height() == evap.height());
  assert(topo.width() == open_water_evap.width() && topo.height() == open_water_evap.height());
  assert(topo.width() == precip.width() && topo.height() == precip.height());
  assert(topo.width() == runoff_ratio.width() && topo.height() == runoff_ratio.height());
  assert(topo.width() == slope.width() && topo.height() == slope.height());
  assert(topo.width() == winter_temp.width() && topo.height() == winter_temp.height());
  assert(topo.width() == ksat.width() && topo.height() == ksat.height());
  assert(topo.width() == land_mask.width() && topo.height() == land_mask.height());
  assert(topo.width() == porosity.width() && topo.height() == porosity.height());
  assert(topo.width() == effective_storativity.width() && topo.height() == effective_storativity.height());
  assert(topo.width() == evap.width() && topo.height() == evap.height());
  assert(topo.width() == fdepth.width() && topo.height() == fdepth.height());
  assert(topo.width() == infiltration_array.width() && topo.height() == infiltration_array.height());
  assert(topo.width() == rech.width() && topo.height() == rech.height());
  assert(topo.width() == runoff.width() && topo.height() == runoff.height());
  assert(topo.width() == transmissivity.width() && topo.height() == transmissivity.height());
  assert(topo.width() == wtd.width() && topo.height() == wtd.height());
  assert(topo.width() == wtd_mid.width() && topo.height() == wtd_mid.height());
  assert(topo.width() == wtd_old.width() && topo.height() == wtd_old.height());
  assert(topo.width() == wtd_T.width() && topo.height() == wtd_T.height());
  assert(topo.width() == wtd_T_iteration.width() && topo.height() == wtd_T_iteration.height());
  assert(topo.width() == original_wtd.width() && topo.height() == original_wtd.height());
  assert(topo.width() == scalar_array_x.width() && topo.height() == scalar_array_x.height());
  assert(topo.width() == scalar_array_y.width() && topo.height() == scalar_array_y.height());
  assert(topo.width() == label.width() && topo.height() == label.height());
  assert(topo.width() == final_label.width() && topo.height() == final_label.height());
  assert(topo.width() == flowdirs.width() && topo.height() == flowdirs.height());
  if (!fdepth_end.empty()) {
    assert(topo.width() == evap_start.width() && topo.height() == evap_start.height());
    assert(topo.width() == evap_end.width() && topo.height() == evap_end.height());
    assert(topo.width() == open_water_evap_start.width() && topo.height() == open_water_evap_start.height());
    assert(topo.width() == open_water_evap_end.width() && topo.height() == open_water_evap_end.height());
    assert(topo.width() == precip_start.width() && topo.height() == precip_start.height());
    assert(topo.width() == precip_end.width() && topo.height() == precip_end.height());
    assert(topo.width() == runoff_ratio_start.width() && topo.height() == runoff_ratio_start.height());
    assert(topo.width() == runoff_ratio_end.width() && topo.height() == runoff_ratio_end.height());
    assert(topo.width() == slope_start.width() && topo.height() == slope_start.height());
    assert(topo.width() == slope_end.width() && topo.height() == slope_end.height());
    assert(topo.width() == topo_start.width() && topo.height() == topo_start.height());
    assert(topo.width() == topo_end.width() && topo.height() == topo_end.height());
    assert(topo.width() == winter_temp_start.width() && topo.height() == winter_temp_start.height());
    assert(topo.width() == winter_temp_end.width() && topo.height() == winter_temp_end.height());
    assert(topo.width() == fdepth_start.width() && topo.height() == fdepth_start.height());
    assert(topo.width() == fdepth_end.width() && topo.height() == fdepth_end.height());
  }
  if (!vert_ksat.empty()) {
    assert(topo.width() == vert_ksat.width() && topo.height() == vert_ksat.height());
  }
}
