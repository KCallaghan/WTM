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

  std::string name        = UNINIT_STR;
  std::string start_name  = UNINIT_STR;
  std::string end_name    = UNINIT_STR;
  std::string surfdatadir = UNINIT_STR;
  std::string region      = UNINIT_STR;
  std::string HAD         = UNINIT_STR;
  std::string run_type    = UNINIT_STR;
  std::string time_start  = UNINIT_STR;
  std::string time_end    = UNINIT_STR;
  std::string textfilename    = UNINIT_STR;

  std::string outfilename    = UNINIT_STR;


 //const int    dltxy = 120; //There are 120 30 arc-second pieces in one degree
  int cells_per_degree = -1;//120;

 // const double dx    = dy;
  const double UNDEF  = -1.0e7;


  bool interpolated = false;
  
  double southern_edge  = std::numeric_limits<double>::signaling_NaN();
  double deltat         = std::numeric_limits<double>::signaling_NaN();
  double wtdmax         = std::numeric_limits<double>::signaling_NaN();
  double cellsize_n_s_metres         = std::numeric_limits<double>::signaling_NaN();
 double infiltration = 0.1;
  

  int total_cycles         = -1;



  //Set for convenience within the code
  int ncells_x  = -1;
  int ncells_y = -1;

  void print() const;
};

#endif
