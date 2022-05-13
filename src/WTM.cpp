#include "fill_spill_merge.hpp"
#include "irf.hpp"
#include "transient_groundwater.hpp"

#include <fmt/core.h>
#include <richdem/common/timer.hpp>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

namespace dh = richdem::dephier;
namespace rd = richdem;

constexpr double seconds_in_a_year = 31536000.;

std::string get_current_time_and_date_as_str() {
  const auto now       = std::chrono::system_clock::now();
  const auto in_time_t = std::chrono::system_clock::to_time_t(now);

  std::stringstream ss;
  ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

static char help[] = "trying petsc method to solve the problem using Newton";

void initialise(Parameters& params, ArrayPack& arp) {
  std::ofstream textfile(params.textfilename, std::ios_base::app);
  // Text file to save outputs of how much is changing and
  // min and max wtd at various times

  if (params.run_type == "transient") {
    textfile << "Initialise transient" << std::endl;
    InitialiseTransient(params, arp);
    // compute changing cell size and distances between cells as
    // these change with latitude:
    cell_size_area(params, arp);
    textfile << "computed distances, areas, and latitudes" << std::endl;
    // finalise some setup for runoff, labels, etc that is the
    // same for both run types.
    InitialiseBoth(params, arp);
  } else if (params.run_type == "equilibrium") {
    textfile << "Initialise equilibrium" << std::endl;
    InitialiseEquilibrium(params, arp);
    // compute changing cell size and distances between cells as
    // these change with latitude:
    cell_size_area(params, arp);
    textfile << "computed distances, areas, and latitudes" << std::endl;
    // finalise some setup for runoff, labels, etc that is the
    // same for both run types.
    InitialiseBoth(params, arp);
  } else if (params.run_type == "test") {
    textfile << "Initialise test" << std::endl;
    InitialiseTest(params, arp);
    // compute changing cell size and distances between cells as
    // these change with latitude:
    cell_size_area(params, arp);
    textfile << "computed distances, areas, and latitudes" << std::endl;
  } else {
    throw std::runtime_error("That was not a recognised run type! Please choose transient or equilibrium.");
  }

  arp.check();

  // Print column headings to textfile to match data that will be printed after each time step.
  textfile << "Cycles_done Total_wtd_change Change_in_GW_only Change_in_SW_only absolute_value_total_wtd_change "
              "abs_change_in_GW abs_change_in_SW change_in_infiltration total_recharge_added total_loss_to_ocean "
              "sum_of_water_tables "
           << std::endl;
  textfile.close();
}

template <class elev_t>
void update(Parameters& params, ArrayPack& arp, richdem::dephier::DepressionHierarchy<elev_t>& deps) {
  richdem::Timer timer_overall;
  timer_overall.start();

  if (params.run_type == "transient") {
    UpdateTransientArrays(params, arp);
    // linear interpolation of input data from start to end times.
    // with transient runs, we have to redo the depression hierarchy every time,
    // since the topography is changing.
    deps = dh::GetDepressionHierarchy<float, rd::Topology::D8>(
        arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);
  }

  // TODO: How should equilibrium know when to exit?
  if ((params.cycles_done % params.cycles_to_save) == 0) {
    // Save the output every "cycles_to_save" iterations, under a new filename
    // so we can compare how the water table has changed through time.
    arp.wtd.saveGDAL(fmt::format("{}{:09}.tif", params.outfile_prefix, params.cycles_done));
  }

  arp.wtd_old = arp.wtd;  // These are used to see how much change occurs
  arp.wtd_mid = arp.wtd;  // in FSM vs in the groundwater portion.

  //////////////////////
  // Move groundwater //
  //////////////////////

  std::cerr << "Before GW time: " << get_current_time_and_date_as_str() << std::endl;

  richdem::Timer time_groundwater;
  time_groundwater.start();

  // These iterations refer to how many times to repeat the time step within the groundwater
  // portion of code before running FSM. For example, 1 year GW then FSM could also be run as
  // 2x 6 months GW then FSM.
  int iter_count = 0;
  while (iter_count++ < params.maxiter) {
    FanDarcyGroundwater::update(params, arp);
  }

  std::cerr << "t GW time = " << time_groundwater.lap() << std::endl;
  std::cerr << "t After GW time: " << get_current_time_and_date_as_str() << std::endl;

  arp.wtd_mid = arp.wtd;

  ////////////////////////
  // Move surface water //
  ////////////////////////

  if (params.fsm_on) {
    richdem::Timer fsm_timer;
    fsm_timer.start();

    dh::FillSpillMerge(params, deps, arp);

    std::cerr << "t FSM time = " << fsm_timer.lap() << std::endl;
    std::cerr << "t After FSM time: " << get_current_time_and_date_as_str() << std::endl;
  }

  /////////////////////////
  // Set recharge values //
  /////////////////////////

  // Check to see where there is surface water, and adjust how evaporation works
  // at these locations.
  richdem::Timer recharge_timer;
  recharge_timer.start();

  // Evap mode 1: Use the computed open-water evaporation rate
  if (params.evap_mode) {
    std::cout << "p updating the recharge field" << std::endl;
#pragma omp parallel for default(none) shared(arp, params)
    for (unsigned int i = 0; i < arp.topo.size(); i++) {
      if (arp.wtd(i) > 0) {  // if there is surface water present
        arp.rech(i) = (arp.precip(i) - arp.open_water_evap(i)) / seconds_in_a_year * params.deltat;
      } else {  // water table is below the surface
        // Recharge is always positive.
        arp.rech(i) =
            (std::max(0., static_cast<double>(arp.precip(i)) - arp.evap(i))) / seconds_in_a_year * params.deltat;
      }

      if (arp.rech(i) > 0) {
        // if there is positive recharge, some of it may run off.
        // set the amount of runoff based on runoff_ratio, and subtract this amount from the recharge.
        arp.runoff(i) = arp.runoff_ratio(i) * arp.rech(i);
        arp.rech(i) -= arp.runoff(i);
      }
    }
  }

  // Evap mode 0: remove all surface water (like Fan Reinfelder et al., 2013)
  else {
    std::cout << "p removing all surface water" << std::endl;
#pragma omp parallel for default(none) shared(arp, params)
    for (unsigned int i = 0; i < arp.topo.size(); i++) {
      if (arp.wtd(i) > 0) {  // if there is surface water present
        arp.wtd(i) = 0;      // use this option when testing GW component alone
        // still set recharge because it could be positive in this cell, and some may run off or move to neighbouring
        // cells
        arp.rech(i) = (arp.precip(i) - arp.open_water_evap(i)) / seconds_in_a_year * params.deltat;
      } else {  // water table is below the surface
        arp.rech(i) =
            (std::max(0., static_cast<double>(arp.precip(i)) - arp.evap(i))) / seconds_in_a_year * params.deltat;
      }
      if (arp.rech(i) > 0) {
        // if there is positive recharge, some of it may run off.
        // set the amount of runoff based on runoff_ratio, and subtract this amount from the recharge.
        arp.runoff(i) = arp.runoff_ratio(i) * arp.rech(i);
        arp.rech(i) -= arp.runoff(i);
      }
    }
  }

  std::cerr << "t Set recharge time = " << recharge_timer.lap() << std::endl;
  std::cerr << "After setting recharge values: " << get_current_time_and_date_as_str() << std::endl;

  // Print values about the change in water table depth to the text file.
  PrintValues(params, arp);

  arp.wtd_old = arp.wtd;
  params.cycles_done += 1;
  std::cerr << "t Done time = " << get_current_time_and_date_as_str() << std::endl;
  std::cerr << "t WTM update time = " << timer_overall.lap() << std::endl;
}

void run(Parameters& params, ArrayPack& arp) {
  // Set the initial depression hierarchy.
  // For equilibrium runs, this is the only time this needs to be done.
  auto deps = dh::GetDepressionHierarchy<float, rd::Topology::D8>(
      arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

  while (params.cycles_done < params.total_cycles) {
    update(params, arp, deps);
  }
}

void finalise(Parameters& params, ArrayPack& arp) {
  std::ofstream textfile(params.textfilename, std::ios_base::app);

  textfile << "p done with processing" << std::endl;
  // save the final answer for water table depth.
  arp.wtd.saveGDAL(fmt::format("{}{:09}.tif", params.outfile_prefix, params.cycles_done));

  textfile.close();
}

int main(int argc, char** argv) {
  if (argc != 2) {
    // Make sure that the user is running the code with a configuration file.
    std::cerr << "Syntax: " << argv[0] << " <Configuration File>" << std::endl;
    return -1;
  }

  std::cerr << "Reading configuration file '" << argv << "'..." << std::endl;
  Parameters params(argv[1]);

  ArrayPack arp;

  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv, (char*)0, help);
  if (ierr)
    return ierr;

  initialise(params, arp);
  run(params, arp);
  finalise(params, arp);

  PetscFinalize();

  return 0;
}
