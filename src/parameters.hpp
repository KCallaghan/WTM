#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <stdint.h>
#include <cmath>
#include <limits>
#include <string>

struct Parameters {
  Parameters() = default;
  Parameters(const std::string& config_file);
  void check() const;
  std::string get_path(const std::string& time, const std::string& layer_name) const;
  std::string get_path(const std::string& layer_name) const;

  static constexpr auto UNINIT_STR = "uninitialized";

  int32_t maxiter = -1;

  std::string outfile_prefix = UNINIT_STR;
  std::string region         = UNINIT_STR;
  std::string run_type       = UNINIT_STR;
  std::string surfdatadir    = UNINIT_STR;
  std::string textfilename   = UNINIT_STR;
  std::string time_start     = UNINIT_STR;
  std::string time_end       = UNINIT_STR;

  double cells_per_degree = -1;

  double UNDEF = -1.0e7;

  int32_t infiltration_on = -1;
  int32_t supplied_wt     = -1;
  int32_t evap_mode       = -1;
  int32_t fsm_on          = -1;
  int32_t runoff_ratio_on = -1;

  double deltat          = std::numeric_limits<double>::signaling_NaN();
  double fdepth_a        = -1.;
  double fdepth_b        = -1.;
  double fdepth_fmin     = -1.;
  double southern_edge   = std::numeric_limits<double>::signaling_NaN();
  int32_t total_cycles   = -1;
  int32_t cycles_to_save = -1;

  double cellsize_n_s_metres = std::numeric_limits<double>::signaling_NaN();
  int32_t cycles_done        = 0;
  double infiltration_change = 0.;

  // Set for convenience within the code
  int32_t ncells_x = -1;
  int32_t ncells_y = -1;

  void print() const;
};

#endif
