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

  int32_t iterations = -1;
  int32_t maxiter    = -1;

  std::string outfile_prefix = UNINIT_STR;
  std::string region         = UNINIT_STR;
  std::string run_type       = UNINIT_STR;
  std::string surfdatadir    = UNINIT_STR;
  std::string textfilename   = UNINIT_STR;
  std::string time_start     = UNINIT_STR;
  std::string time_end       = UNINIT_STR;

  double cells_per_degree = -1;

  double UNDEF = -1.0e7;

  bool infiltration_on;
  bool supplied_wt;
  bool evap_mode = 1;  // Default potential evaporation
  bool fsm_on    = 1;  // Default surface water on

  double deltat             = std::numeric_limits<double>::signaling_NaN();
  double fdepth_a           = 0.;
  double fdepth_b           = 0.;
  double fdepth_fmin        = 0.;
  double southern_edge      = std::numeric_limits<double>::signaling_NaN();
  int32_t picard_iterations = 1;  // Default only one iteration
  int32_t total_cycles      = -1;
  int32_t parallel_threads  = 1;  // one thread by default to accommodate computer architecture

  double cellsize_n_s_metres  = std::numeric_limits<double>::signaling_NaN();
  int32_t cycles_done         = 0;
  double infiltration_change  = 0.;
  double x_partial            = 0.;
  double total_added_recharge = 0.;
  double total_loss_to_ocean  = 0.;

  // Set for convenience within the code
  int32_t ncells_x = -1;
  int32_t ncells_y = -1;

  void print() const;
};

#endif
