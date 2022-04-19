#include "parameters.hpp"

#include <fmt/core.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Real initializer
Parameters::Parameters(const std::string& config_file) {
  std::ifstream fin(config_file);

  if (!fin.good()) {
    throw std::runtime_error("Failed to read config file!");
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) {
      continue;
    }

    std::stringstream ss(line);
    std::string key;
    ss >> key;

    // Dummy key to make it easier to alphabetize list below
    if (key.empty()) {
    } else if (key == "cells_per_degree") {
      ss >> cells_per_degree;
    } else if (key == "cycles_to_save") {
      ss >> cycles_to_save;
    } else if (key == "deltat") {
      ss >> deltat;
    } else if (key == "evap_mode") {
      ss >> evap_mode;
    } else if (key == "fdepth_a") {
      ss >> fdepth_a;
    } else if (key == "fdepth_b") {
      ss >> fdepth_b;
    } else if (key == "fdepth_fmin") {
      ss >> fdepth_fmin;
    } else if (key == "fsm_on") {
      ss >> fsm_on;
    } else if (key == "infiltration_on") {
      ss >> infiltration_on;
    } else if (key == "maxiter") {
      ss >> maxiter;
    } else if (key == "outfile_prefix") {
      ss >> outfile_prefix;
    } else if (key == "parallel_threads") {
      ss >> parallel_threads;
    } else if (key == "picard_iterations") {
      ss >> picard_iterations;
    } else if (key == "region") {
      ss >> region;
    } else if (key == "run_type") {
      ss >> run_type;
    } else if (key == "runoff_ratio_on") {
      ss >> runoff_ratio_on;
    } else if (key == "solver_tolerance_value") {
      ss >> solver_tolerance_value;
    } else if (key == "southern_edge") {
      ss >> southern_edge;
    } else if (key == "supplied_wt") {
      ss >> supplied_wt;
    } else if (key == "surfdatadir") {
      ss >> surfdatadir;
    } else if (key == "textfilename") {
      ss >> textfilename;
    } else if (key == "time_end") {
      ss >> time_end;
    } else if (key == "time_start") {
      ss >> time_start;
    } else if (key == "total_cycles") {
      ss >> total_cycles;
    } else {
      throw std::runtime_error("Unrecognised key: " + key);
    }
  }
  std::cout << infiltration_on << std::endl;

  check();
}

void Parameters::check() const {
  if (cells_per_degree < 0) {
    throw std::runtime_error("please enter a positive value for cells_per_degree!");
  }
  if (cycles_to_save < 0) {
    throw std::runtime_error("please enter a positive value for cycles_to_save!");
  }
  if (std::isnan(deltat) || deltat < 0) {
    throw std::runtime_error("please enter a positive value for deltat!");
  }
  if (std::isnan(southern_edge) || southern_edge < -90 || southern_edge > 90) {
    throw std::runtime_error("please enter a value between -90 and 90 degrees for the southern_edge!");
  }
  if (evap_mode != 0 && evap_mode != 1) {
    throw std::runtime_error(
        "set evap_mode to 0 to remove all surface water, or 1 to use a grid of potential evaporation for lakes.");
  }
  if (fdepth_a < 0) {
    throw std::runtime_error("please enter a positive value for fdepth_a!");
  }
  if (fdepth_b < 0) {
    throw std::runtime_error("please enter a positive value for fdepth_b!");
  }
  if (fdepth_fmin < 0) {
    throw std::runtime_error("please enter a positive value for fdepth_fmin!");
  }
  if (picard_iterations < 0) {
    throw std::runtime_error("please enter a positive value for picard_iterations!");
  }
  if (fsm_on != 0 && fsm_on != 1) {
    throw std::runtime_error(
        "set fsm_on to 1 to allow Fill-Spill-Merge to move surface water, or 0 to disable Fill-Spill-Merge.");
  }
  if (infiltration_on != 0 && infiltration_on != 1) {
    throw std::runtime_error(
        "set infiltration_on to 1 to allow water to infiltrate as it flows downslope, or 0 to neglect infiltration and "
        "assume impermeable substrates while flowing downslope.");
  }
  if (runoff_ratio_on != 0 && runoff_ratio_on != 1) {
    throw std::runtime_error(
        "set runoff_ratio_on to 1 to supply a runoff ratio array, or 0 to assume all P-ET infiltrates in the cell "
        "where it falls.");
  }
  if (solver_tolerance_value < 0) {
    throw std::runtime_error("select a positive value for solver_tolerance_value");
  }
  if (supplied_wt != 0 && supplied_wt != 1) {
    throw std::runtime_error(
        "set supplied_wt to 1 to supply a starting water table, or 0 to set starting water table == 0 (only available "
        "for equilibrium runs).");
  }
  if (maxiter == -1 || maxiter != static_cast<int32_t>(maxiter)) {
    throw std::runtime_error("please enter a positive integer value for maxiter!");
  }
  if (outfile_prefix == UNINIT_STR) {
    throw std::runtime_error("please provide an outfile_prefix!");
  }
  if (region == UNINIT_STR) {
    throw std::runtime_error("please provide a region!");
  }
  if (run_type == UNINIT_STR) {
    throw std::runtime_error("please provide a run_type!");
  }
  if (surfdatadir == UNINIT_STR) {
    throw std::runtime_error("please provide a surfdatadir!");
  }
  if (textfilename == UNINIT_STR) {
    throw std::runtime_error("please provide a textfilename!");
  }
  if (time_start == UNINIT_STR) {
    throw std::runtime_error("please provide a time_start!");
  }
  if (time_end == UNINIT_STR) {
    throw std::runtime_error("please provide a time_end!");
  }
  if (parallel_threads == -1) {
    throw std::runtime_error("please enter a positive integer value for parallel_threads!");
  }
  if (total_cycles == -1) {
    throw std::runtime_error("please enter a positive integer value for total_cycles!");
  }
}

std::string Parameters::get_path(const std::string& time, const std::string& layer_name) const {
  constexpr auto SURF_DATA_PATH_FORMAT = "{}/{}_{}_{}.tif";
  return fmt::format(SURF_DATA_PATH_FORMAT, surfdatadir, region, time, layer_name);
}

std::string Parameters::get_path(const std::string& layer_name) const {
  constexpr auto SURF_DATA_PATH_FORMAT = "{}/{}_{}.tif";
  return fmt::format(SURF_DATA_PATH_FORMAT, surfdatadir, region, layer_name);
}

void Parameters::print() const {
  std::cout << "c cells_per_degree       = " << cells_per_degree << std::endl;
  std::cout << "c cycles_to_save         = " << cycles_to_save << std::endl;
  std::cout << "c deltat                 = " << deltat << std::endl;
  std::cout << "c evap_mode              = " << evap_mode << std::endl;
  std::cout << "c fdepth_a               = " << fdepth_a << std::endl;
  std::cout << "c fdepth_b               = " << fdepth_b << std::endl;
  std::cout << "c fdepth_fmin            = " << fdepth_fmin << std::endl;
  std::cout << "c fsm_on                 = " << fsm_on << std::endl;
  std::cout << "c infiltration_on        = " << infiltration_on << std::endl;
  std::cout << "c maxiter                = " << maxiter << std::endl;
  std::cout << "c outfile_prefix         = " << outfile_prefix << std::endl;
  std::cout << "c parallel_threads       = " << parallel_threads << std::endl;
  std::cout << "c picard_iterations      = " << picard_iterations << std::endl;
  std::cout << "c region                 = " << region << std::endl;
  std::cout << "c run_type               = " << run_type << std::endl;
  std::cout << "c runoff_ratio_on        = " << runoff_ratio_on << std::endl;
  std::cout << "c southern_edge          = " << southern_edge << std::endl;
  std::cout << "c solver_tolerance_value = " << solver_tolerance_value << std::endl;
  std::cout << "c supplied_wt            = " << supplied_wt << std::endl;
  std::cout << "c surfdatadir            = " << surfdatadir << std::endl;
  std::cout << "c textfilename           = " << textfilename << std::endl;
  std::cout << "c time_end               = " << time_end << std::endl;
  std::cout << "c time_start             = " << time_start << std::endl;
  std::cout << "c total_cycles           = " << total_cycles << std::endl;
}
