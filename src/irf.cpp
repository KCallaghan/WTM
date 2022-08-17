#include "ArrayPack.hpp"
#include "fill_spill_merge.hpp"
#include "parameters.hpp"

#include <richdem/common/Array2D.hpp>

namespace rd = richdem;
namespace dh = richdem::dephier;

constexpr double UNDEF             = -1.0e7;
constexpr double seconds_in_a_year = 31536000.;

/// We calculate the e-folding depth here, using temperature and slope.
double setup_fdepth(const Parameters& params, const double slope, const double temperature) {
  const auto fdepth = std::max(params.fdepth_a / (1 + params.fdepth_b * slope), params.fdepth_fmin);
  if (temperature > -5) {  // then fdepth = f from Ying's equation S7.
    return fdepth;
  } else if (temperature < -14) {  // then fdpth = f*fT, Ying's equations S7 and S8.
    return fdepth * std::max(0.05, 0.17 + 0.005 * temperature);
  } else {
    return fdepth * std::min(1.0, 1.5 + 0.1 * temperature);
  }
}

/// This function initialises those arrays that are needed only for transient
/// model runs. This includes both start and end states for slope, precipitation,
/// temperature, topography, ET, and relative humidity. We also have a land vs
/// ocean mask for the end time. It also includes the starting water table depth
/// array, a requirement for transient runs.
void InitialiseTransient(Parameters& params, ArrayPack& arp) {
  // width and height in number of cells in the array

  arp.topo_start = rd::Array2D<float>(params.get_path(params.time_start, "topography"));

  params.ncells_x = arp.topo_start.width();
  params.ncells_y = arp.topo_start.height();

  arp.slope_start           = rd::Array2D<float>(params.get_path(params.time_start, "slope"));
  arp.precip_start          = rd::Array2D<float>(params.get_path(params.time_start, "precipitation"));
  arp.evap_start            = rd::Array2D<float>(params.get_path(params.time_start, "evaporation"));
  arp.open_water_evap_start = rd::Array2D<float>(params.get_path(params.time_start, "open_water_evaporation"));
  arp.winter_temp_start     = rd::Array2D<float>(params.get_path(params.time_start, "winter_temperature"));
  arp.topo_end              = rd::Array2D<float>(params.get_path(params.time_end, "topography"));
  arp.slope_end             = rd::Array2D<float>(params.get_path(params.time_end, "slope"));
  arp.land_mask             = rd::Array2D<float>(params.get_path(params.time_end, "mask"));

  // there is land and 0 in the ocean
  arp.land_mask.setEdges(0);

  arp.precip_end          = rd::Array2D<float>(params.get_path(params.time_end, "precipitation"));
  arp.evap_end            = rd::Array2D<float>(params.get_path(params.time_end, "evaporation"));
  arp.open_water_evap_end = rd::Array2D<float>(params.get_path(params.time_end, "open_water_evaporation"));
  arp.winter_temp_end     = rd::Array2D<float>(params.get_path(params.time_end, "winter_temperature"));

  if (params.runoff_ratio_on) {
    arp.runoff_ratio_start = rd::Array2D<float>(params.get_path(params.time_start, "runoff_ratio"));
    arp.runoff_ratio_end   = rd::Array2D<float>(params.get_path(params.time_end, "runoff_ratio"));
  } else {
    arp.runoff_ratio_start = rd::Array2D<float>(arp.topo, 0.0);
    arp.runoff_ratio_end   = rd::Array2D<float>(arp.topo, 0.0);
  }

  if (params.infiltration_on) {
    arp.vert_ksat = rd::Array2D<float>(params.get_path("vertical_ksat"));
  }

  arp.scalar_array_x = rd::Array2D<double>(arp.topo, 0.0);
  arp.scalar_array_y = rd::Array2D<double>(arp.topo, 0.0);

  // load in the wtd result from the previous time:
  arp.wtd = rd::Array2D<double>(params.get_path(params.time_start, "wtd"));

  // calculate the fdepth (e-folding depth, representing rate of decay of the
  // hydraulic conductivity with depth) arrays:
  arp.fdepth_start = rd::Array2D<double>(arp.topo_start, 0);
  arp.fdepth_end   = rd::Array2D<double>(arp.topo_start, 0);

  for (size_t i = 0; i < arp.topo_start.size(); i++) {
    arp.fdepth_start(i) = setup_fdepth(params, arp.slope_start(i), arp.winter_temp_start(i));
    arp.fdepth_end(i)   = setup_fdepth(params, arp.slope_end(i), arp.winter_temp_end(i));
  }

  // initialise the arrays to be as at the starting time:
  arp.topo            = arp.topo_start;
  arp.slope           = arp.slope_start;
  arp.precip          = arp.precip_start;
  arp.runoff_ratio    = arp.runoff_ratio_start;
  arp.evap            = arp.evap_start;
  arp.open_water_evap = arp.open_water_evap_start;
  arp.winter_temp     = arp.winter_temp_start;
  arp.fdepth          = arp.fdepth_start;
  arp.wtd_T           = arp.wtd;
  arp.wtd_T_iteration = arp.wtd;
  arp.original_wtd    = arp.wtd;
  arp.scalar_array_x  = rd::Array2D<double>(arp.topo, 0.0);
  arp.scalar_array_y  = rd::Array2D<double>(arp.topo, 0.0);
}

/// This function initialises those arrays that are needed only for equilibrium
/// model runs.
/// This includes a single array for each of slope, precipitation, temperature,
/// topography, ET, land vs ocean mask, and relative humidity.
/// It also includes setting the starting water table depth array to
/// zero everywhere.
void InitialiseEquilibrium(Parameters& params, ArrayPack& arp) {
  arp.topo = rd::Array2D<float>(params.get_path(params.time_start, "topography"));

  // width and height in number of cells in the array
  params.ncells_x = arp.topo.width();
  params.ncells_y = arp.topo.height();

  arp.slope     = rd::Array2D<float>(params.get_path(params.time_start, "slope"));
  arp.land_mask = rd::Array2D<float>(
      params.get_path(params.time_start, "mask"));  // A binary mask that is 1 where there is land and 0 in the ocean
  arp.land_mask.setEdges(0);

  arp.precip = rd::Array2D<float>(params.get_path(params.time_start, "precipitation"));  // Units: m/yr.
  arp.evap   = rd::Array2D<float>(params.get_path(params.time_start, "evaporation"));    // Units: m/yr.
  arp.open_water_evap =
      rd::Array2D<float>(params.get_path(params.time_start, "open_water_evaporation"));  // Units: m/yr.
  arp.winter_temp =
      rd::Array2D<float>(params.get_path(params.time_start, "winter_temperature"));  // Units: degrees Celsius

  arp.scalar_array_x = rd::Array2D<double>(arp.topo, 0.0);
  arp.scalar_array_y = rd::Array2D<double>(arp.topo, 0.0);

  if (params.runoff_ratio_on) {
    arp.runoff_ratio = rd::Array2D<float>(params.get_path(params.time_start, "runoff_ratio"));  // Units: m/yr.
  } else {
    arp.runoff_ratio = rd::Array2D<float>(arp.topo, 0.0);  // Units: m/yr.
  }

  if (params.infiltration_on == true) {
    arp.vert_ksat = rd::Array2D<float>(params.get_path("vertical_ksat"));  // Units of ksat are m/s.
  }

  if (params.supplied_wt == true) {
    arp.wtd             = rd::Array2D<double>(params.get_path(params.time_start, "starting_wt"));
    arp.wtd_T           = arp.wtd;
    arp.wtd_T_iteration = arp.wtd;
    arp.original_wtd    = arp.wtd;
  } else {
    arp.wtd             = rd::Array2D<double>(arp.topo, 0.);
    arp.wtd_T           = rd::Array2D<double>(arp.topo, 0.);
    arp.wtd_T_iteration = rd::Array2D<double>(arp.topo, 0.);
    arp.original_wtd    = rd::Array2D<double>(arp.topo, 0.);
  }
  // we start with a water table at the surface for equilibrium runs.

  arp.fdepth = rd::Array2D<double>(arp.topo, 0);
  for (unsigned int i = 0; i < arp.topo.size(); i++) {
    arp.fdepth(i) = setup_fdepth(params, arp.slope(i), arp.winter_temp(i));
  }
}

void InitialiseTest(Parameters& params, ArrayPack& arp) {
  arp.topo  = rd::Array2D<float>(params.get_path("topography"));
  arp.slope = rd::Array2D<float>(params.get_path("slope"));  // Slope as a value from 0 to 1.

  // width and height in number of cells in the array
  params.ncells_x = arp.topo.width();
  params.ncells_y = arp.topo.height();

  if (params.infiltration_on) {
    arp.vert_ksat = rd::Array2D<float>(arp.topo, 0.00001f);  // Units of ksat are m/s.
  }

  // A binary mask that is 1 where there is land and 0 in the ocean
  arp.land_mask = rd::Array2D<uint8_t>(arp.topo, 1);
  arp.land_mask.setEdges(0);

  arp.precip          = rd::Array2D<float>(arp.topo, 0.03f);  // Units: m/yr.
  arp.runoff_ratio    = rd::Array2D<float>(arp.topo, 0.);     // Units: m/yr.
  arp.evap            = rd::Array2D<float>(arp.topo, 0.);     // Units: m/yr.
  arp.open_water_evap = rd::Array2D<float>(arp.topo, 0.5);    // Units: m/yr.

  arp.winter_temp     = rd::Array2D<float>(arp.topo, 0);  // Units: deg C
  arp.wtd             = rd::Array2D<double>(arp.topo, 0.0);
  arp.wtd_T           = rd::Array2D<double>(arp.topo, 0.0);
  arp.wtd_T_iteration = rd::Array2D<double>(arp.topo, 0.0);
  arp.original_wtd    = rd::Array2D<double>(arp.topo, 0.0);

  arp.scalar_array_x = rd::Array2D<double>(arp.topo, 0.0);
  arp.scalar_array_y = rd::Array2D<double>(arp.topo, 0.0);

  arp.fdepth = rd::Array2D<double>(arp.topo, 60.);

  // border of 'ocean' with land everywhere else
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask.isEdgeCell(x, y)) {
        arp.land_mask(x, y) = 0;
      } else {
        arp.land_mask(x, y) = 1;
        if (std::isnan(arp.topo(x, y))) {
          arp.topo(x, y) = 0;
        }
      }
    }
  }

  arp.ksat                  = rd::Array2D<float>(arp.topo, 0.0001f);  // Units of ksat are m/s.
  arp.porosity              = rd::Array2D<float>(arp.topo, 0.25);     // Units: unitless
  arp.effective_storativity = rd::Array2D<double>(arp.topo, 0.25);

  // Set arrays that start off with zero or other values,
  // that are not imported files. Just to initialise these -
  // we'll add the appropriate values later.

  // These two are just informational, to see how much change
  // happens in FSM vs in groundwater
  arp.wtd_old = arp.wtd;
  arp.wtd_mid = arp.wtd;

  arp.runoff = rd::Array2D<double>(arp.ksat, 0);

  // This is used to see how much change occurred in infiltration
  // portion of the code. Just informational.
  arp.infiltration_array = rd::Array2D<double>(arp.ksat, 0);

  arp.rech           = rd::Array2D<double>(arp.ksat, 0);
  arp.transmissivity = rd::Array2D<double>(arp.ksat, 0);

  // These are populated during the calculation of the depression hierarchy:
  // No cells are part of a depression
  arp.label = rd::Array2D<dh::dh_label_t>(params.ncells_x, params.ncells_y, dh::NO_DEP);
  // No cells are part of a depression
  arp.final_label = rd::Array2D<dh::dh_label_t>(params.ncells_x, params.ncells_y, dh::NO_DEP);
  // No cells flow anywhere
  arp.flowdirs = rd::Array2D<rd::flowdir_t>(params.ncells_x, params.ncells_y, rd::NO_FLOW);

  // Change undefined cells to 0
  for (size_t i = 0; i < arp.topo.size(); i++) {
    if (arp.topo(i) <= UNDEF) {
      arp.topo(i) = 0;
    }
  }

// get the starting runoff using precip and evap inputs:
#pragma omp parallel for default(none) shared(arp, params)
  for (size_t i = 0; i < arp.topo.size(); i++) {
    arp.rech(i) = (std::max(0., static_cast<double>(arp.precip(i)) - arp.evap(i))) / seconds_in_a_year * params.deltat;
    if (arp.porosity(i) <= 0) {
      arp.porosity(i) = 0.0000001f;  // not sure why it is sometimes processing cells with 0 porosity?
    }
  }

// Wtd is 0 in the ocean and under the ice:
#pragma omp parallel for default(none) shared(arp)
  for (size_t i = 0; i < arp.topo.size(); i++) {
    if (arp.land_mask(i) == 0) {  // || arp.ice_mask(i) == 1){
      arp.wtd(i)             = 0.;
      arp.wtd_T(i)           = 0.;
      arp.wtd_T_iteration(i) = 0.;
      arp.original_wtd(i)    = 0.;
      // arp.topo(i)            = 0.;
    }
  }

// Label the ocean cells. This is a precondition for
// using `GetDepressionHierarchy()`.
#pragma omp parallel for default(none) shared(arp)
  for (unsigned int i = 0; i < arp.label.size(); i++) {
    if (arp.land_mask(i) == 0) {
      arp.label(i)       = dh::OCEAN;
      arp.final_label(i) = dh::OCEAN;
    }
  }
}

/// Here, we use the number of cells per degree (a user-defined value),
/// the southern-most latitude of the domain (also user-defined),
/// and the radius of the Earth to calculate the latitude of each row of cells,
/// the size of a cell in the N-S and E-W directions, and the area of each cell.
void cell_size_area(Parameters& params, ArrayPack& arp) {
  // compute changing cell size and distances between
  // cells as these change with latitude:

  constexpr double earth_radius = 6371000.;    // metres
  constexpr double deg_to_rad   = M_PI / 180;  // convert degrees to radians

  // Radius * Pi = Distance from N to S pole
  // Distance / 180 = Meters / degree latitude
  const auto meters_per_degree = earth_radius * deg_to_rad;

  // N-S Meters per cell; distance between lines of latitude is a constant
  params.cellsize_n_s_metres = meters_per_degree / params.cells_per_degree;

  // initialise some arrays

  // size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cellsize_e_w_metres.resize(params.ncells_y);
  // cell area (metres squared)
  arp.cell_area.resize(params.ncells_y);

  // used to calculate cell latitude in radians.
  // southern edge of the domain in degrees, plus the number of cells up from this
  // location/the number of cells per degree, converted to radians.
  const auto cell_position_latitude = [&](const auto cell_idx) {
    return (cell_idx / params.cells_per_degree + params.southern_edge) * deg_to_rad;
  };

  for (int j = 0; j < params.ncells_y; j++) {
    // latitude, in radians, at the southern edge of a cell:
    const double latitude_radians_S = cell_position_latitude(j);
    // latitude, in radians, at the northern edge of a cell (add a cell; equal to the southern edge of the next cell):
    const double latitude_radians_N = cell_position_latitude(j + 1);

    // distance between lines of longitude varies with latitude.

    // distance at the northern edge of the cell for the given latitude:
    const double cellsize_e_w_metres_N = params.cellsize_n_s_metres * std::cos(latitude_radians_N);
    // distance at the southern edge of the cell for the given latitude:
    const double cellsize_e_w_metres_S = params.cellsize_n_s_metres * std::cos(latitude_radians_S);

    arp.cellsize_e_w_metres[j] = (cellsize_e_w_metres_N + cellsize_e_w_metres_S) / 2.;

    // cell area computed as a trapezoid, using unchanging north-south distance,
    // and east-west average distance.
    arp.cell_area[j] = params.cellsize_n_s_metres * arp.cellsize_e_w_metres[j];

    if (arp.cell_area[j] < 0) {
      throw std::runtime_error("Cell with a negative area was found!");
    }
  }
}

/// This function initialises those arrays that are used for both equilibrium
/// and transient model runs. This includes arrays that start off with zero
/// values, as well as the label, final_label, and flowdirs arrays.
void InitialiseBoth(const Parameters& params, ArrayPack& arp) {
  arp.ksat     = rd::Array2D<float>(params.get_path("horizontal_ksat"));
  arp.porosity = rd::Array2D<float>(params.get_path("porosity"));

  arp.effective_storativity = rd::Array2D<double>(arp.topo, 0.);
  // Set arrays that start off with zero or other values,
  // that are not imported files. Just to initialise these -
  // we'll add the appropriate values later.

  // These two are just informational, to see how much change
  // happens in FSM vs in groundwater
  arp.wtd_old = arp.wtd;
  arp.wtd_mid = arp.wtd;

  arp.runoff = rd::Array2D<double>(arp.ksat, 0);

  // These are used to see how much change occurred in infiltration
  // and updating lakes portions of the code. Just informational.
  arp.infiltration_array = rd::Array2D<double>(arp.ksat, 0);

  arp.rech           = rd::Array2D<double>(arp.ksat, 0);
  arp.transmissivity = rd::Array2D<double>(arp.ksat, 0);

  // These are populated during the calculation of the depression hierarchy:
  // No cells are part of a depression
  arp.label = rd::Array2D<dh::dh_label_t>(params.ncells_x, params.ncells_y, dh::NO_DEP);
  // No cells are part of a depression
  arp.final_label = rd::Array2D<dh::dh_label_t>(params.ncells_x, params.ncells_y, dh::NO_DEP);
  // No cells flow anywhere
  arp.flowdirs = rd::Array2D<rd::flowdir_t>(params.ncells_x, params.ncells_y, rd::NO_FLOW);

  // Change undefined cells to 0
  for (unsigned int i = 0; i < arp.topo.size(); i++) {
    if (arp.topo(i) <= UNDEF) {
      arp.topo(i) = 0;
    }
  }

  // get the starting runoff using precip and evap inputs:
#pragma omp parallel for default(none) shared(arp)
  for (unsigned int i = 0; i < arp.topo.size(); i++) {
    arp.rech(i) = (std::max(0., static_cast<double>(arp.precip(i)) - arp.evap(i))) / seconds_in_a_year * params.deltat;
    if (arp.porosity(i) <= 0) {
      arp.porosity(i) = 0.0000001f;  // not sure why it is sometimes processing cells with 0 porosity?
    }
  }

  // Wtd is 0 in the ocean and under the ice:
#pragma omp parallel for default(none) shared(arp)
  for (unsigned int i = 0; i < arp.topo.size(); i++) {
    if (arp.land_mask(i) == 0) {  //|| arp.ice_mask(i) ==1){
      arp.wtd(i)             = 0.;
      arp.wtd_T(i)           = 0.;
      arp.wtd_T_iteration(i) = 0.;
      arp.original_wtd(i)    = 0.;
      arp.topo(i)            = 0.;
    }
  }

// Label the ocean cells. This is a precondition for
// using `GetDepressionHierarchy()`.
#pragma omp parallel for default(none) shared(arp)
  for (unsigned int i = 0; i < arp.label.size(); i++) {
    if (arp.land_mask(i) == 0) {
      arp.label(i)       = dh::OCEAN;
      arp.final_label(i) = dh::OCEAN;
    }
  }
}

/// In transient runs, we adjust the input arrays via a
// linear interpolation from the start state to the end state at each iteration.
/// We do so here, and also reset the label and flow direction arrays,
/// since the depression hierarchy needs to be
/// recalculated due to the changed topography.
void UpdateTransientArrays(const Parameters& params, ArrayPack& arp) {
  for (unsigned int i = 0; i < arp.topo.size(); i++) {
    const auto f = params.cycles_done / params.total_cycles;

    arp.topo(i)            = (1 - f) * arp.topo_start(i) + f * arp.topo_end(i);
    arp.slope(i)           = (1 - f) * arp.slope_start(i) + f * arp.slope_end(i);
    arp.precip(i)          = (1 - f) * arp.precip_start(i) + f * arp.precip_end(i);
    arp.runoff_ratio(i)    = (1 - f) * arp.runoff_ratio_start(i) + f * arp.runoff_ratio_end(i);
    arp.evap(i)            = (1 - f) * arp.evap_start(i) + f * arp.evap_end(i);
    arp.open_water_evap(i) = (1 - f) * arp.open_water_evap_start(i) + f * arp.open_water_evap_end(i);
    arp.winter_temp(i)     = (1 - f) * arp.winter_temp_start(i) + f * arp.winter_temp_end(i);
    arp.fdepth(i)          = (1 - f) * arp.fdepth_start(i) + f * arp.fdepth_end(i);

    arp.label(i)       = dh::NO_DEP;   // No cells are part of a depression
    arp.final_label(i) = dh::NO_DEP;   // No cells are part of a depression
    arp.flowdirs(i)    = rd::NO_FLOW;  // No cells flow anywhere
  }

#pragma omp parallel for default(none) shared(arp)
  for (unsigned int i = 0; i < arp.label.size(); i++) {
    if (arp.land_mask(i) == 0) {
      arp.label(i)       = dh::OCEAN;
      arp.final_label(i) = dh::OCEAN;
    }
  }
}

/// In this function, we use a few of the variables that were created for
/// informational purposes to help us understand how much the water table
/// is changing per iteration, and where in
/// the code that change is occurring. We print these values to a text file.
void PrintValues(Parameters& params, const ArrayPack& arp) {
  std::ofstream textfile(params.textfilename, std::ios_base::app);

  double abs_total_wtd_change = 0.0;
  double abs_wtd_mid_change   = 0.0;
  double abs_GW_wtd_change    = 0.0;
  double total_wtd_change     = 0.0;
  double wtd_mid_change       = 0.0;
  double GW_wtd_change        = 0.0;
  double wtd_sum              = 0.0;

  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      abs_total_wtd_change += std::abs(arp.wtd(x, y) - arp.wtd_old(x, y));
      abs_wtd_mid_change += std::abs(arp.wtd(x, y) - arp.wtd_mid(x, y));
      abs_GW_wtd_change += std::abs(arp.wtd_mid(x, y) - arp.wtd_old(x, y));
      total_wtd_change += (arp.wtd(x, y) - arp.wtd_old(x, y));
      wtd_mid_change += (arp.wtd(x, y) - arp.wtd_mid(x, y));
      GW_wtd_change += (arp.wtd_mid(x, y) - arp.wtd_old(x, y));
      params.infiltration_change += arp.infiltration_array(x, y);
      if (arp.wtd(x, y) > 0) {
        wtd_sum += arp.wtd(x, y) * arp.cell_area[y];
      } else {
        wtd_sum += arp.wtd(x, y) * arp.porosity(x, y) * arp.cell_area[y];
      }
    }
  }

  textfile << params.cycles_done << " " << total_wtd_change << " " << GW_wtd_change << " " << wtd_mid_change << " "
           << abs_total_wtd_change << " " << abs_GW_wtd_change << " " << abs_wtd_mid_change << " "
           << params.infiltration_change << " " << arp.total_added_recharge << " " << arp.total_loss_to_ocean << " "
           << wtd_sum << " " << std::endl;

  textfile.close();
}
