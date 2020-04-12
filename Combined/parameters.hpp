#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <cmath>
#include <limits>
#include <string>

const std::string UNINIT_STR = "uninitialized";

class Parameters {
 public:
  Parameters(const std::string config_file);

  int iterations = -1;
  int maxiter    = -1;

  std::string surfdatadir  = UNINIT_STR;
  std::string region       = UNINIT_STR;
  std::string run_type     = UNINIT_STR;
  std::string time_start   = UNINIT_STR;
  std::string time_end     = UNINIT_STR;
  std::string textfilename = UNINIT_STR;
  std::string outfilename  = UNINIT_STR;

  int cells_per_degree = -1;

  const double UNDEF  = -1.0e7;

  bool infiltration_on;
  
  double southern_edge        = std::numeric_limits<double>::signaling_NaN();
  double deltat               = std::numeric_limits<double>::signaling_NaN();
  double cellsize_n_s_metres  = std::numeric_limits<double>::signaling_NaN();
  float  infiltration         = 0.0;
  int    cycles_done          = 0;
  float  total_wtd_change     = 0.0;
  float  wtd_mid_change       = 0.0;
  float  GW_wtd_change        = 0.0;
  float  abs_total_wtd_change = 0.0;
  float  abs_wtd_mid_change   = 0.0;
  float  abs_GW_wtd_change    = 0.0;
  float  infiltration_change  = 0.0;
  float  surface_change       = 0.0;
  int    total_cycles         = -1;

  //Set for convenience within the code
  int ncells_x  = -1;
  int ncells_y  = -1;

  void print() const;
};

#endif
