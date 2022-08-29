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
    } else if (key == "region") {
      ss >> region;
    } else if (key == "run_type") {
      ss >> run_type;
    } else if (key == "runoff_ratio_on") {
      ss >> runoff_ratio_on;
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
  const auto check_positive = [](const std::string name, const auto val) {
    if (std::isnan(val) || val < 0) {
      throw std::runtime_error("Please enter a positive value for " + name);
    }
  };

  const auto check_string_init = [&](const std::string name, const std::string& val) {
    if (val == UNINIT_STR) {
      throw std::runtime_error("Please provide a value for " + name);
    }
  };

  const auto check_binary = [](const auto val, const std::string& msg) {
    if (val != 0 && val != 1) {
      throw std::runtime_error(msg);
    }
  };

  check_positive("cells_per_degree", cells_per_degree);
  check_positive("cycles_to_save", cycles_to_save);
  check_positive("deltat", deltat);
  if (std::isnan(southern_edge) || southern_edge < -90 || southern_edge > 90) {
    throw std::runtime_error("please enter a value between -90 and 90 degrees for the southern_edge!");
  }
  check_binary(
      evap_mode,
      "set evap_mode to 0 to remove all surface water, or 1 to use a grid of potential evaporation for lakes.");
  check_positive("fdepth_a", fdepth_a);
  check_positive("fdepth_b", fdepth_b);
  check_positive("fdepth_fmin", fdepth_fmin);
  check_binary(
      fsm_on, "set fsm_on to 1 to allow Fill-Spill-Merge to move surface water, or 0 to disable Fill-Spill-Merge.");
  check_binary(
      infiltration_on,
      "set infiltration_on to 1 to allow water to infiltrate as it flows downslope, or 0 to neglect infiltration and "
      "assume impermeable substrates while flowing downslope.");
  check_binary(
      runoff_ratio_on,
      "set runoff_ratio_on to 1 to supply a runoff ratio array, or 0 to assume all P-ET infiltrates in the cell "
      "where it falls.");
  check_binary(
      supplied_wt,
      "set supplied_wt to 1 to supply a starting water table, or 0 to set starting water table == 0 (only available "
      "for equilibrium runs).");
  if (maxiter == -1 || maxiter != static_cast<int32_t>(maxiter)) {
    throw std::runtime_error("please enter a positive integer value for maxiter!");
  }
  check_string_init("outfile_prefix", outfile_prefix);
  check_string_init("region", region);
  check_string_init("run_type", run_type);
  check_string_init("surfdatadir", surfdatadir);
  check_string_init("textfilename", textfilename);
  check_string_init("time_start", time_start);
  check_string_init("time_end", time_end);
  check_positive("total_cycles", total_cycles);
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
  std::cout << "c region                 = " << region << std::endl;
  std::cout << "c run_type               = " << run_type << std::endl;
  std::cout << "c runoff_ratio_on        = " << runoff_ratio_on << std::endl;
  std::cout << "c southern_edge          = " << southern_edge << std::endl;
  std::cout << "c supplied_wt            = " << supplied_wt << std::endl;
  std::cout << "c surfdatadir            = " << surfdatadir << std::endl;
  std::cout << "c textfilename           = " << textfilename << std::endl;
  std::cout << "c time_end               = " << time_end << std::endl;
  std::cout << "c time_start             = " << time_start << std::endl;
  std::cout << "c total_cycles           = " << total_cycles << std::endl;
}
