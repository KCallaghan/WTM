#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <string>
#include <limits>

class Parameters {
 public:
  int iterations = -1;
  int width      = -1;
  int height     = -1;

  std::string name;

  std::string start_name;
  std::string end_name;

  std::string surfdatadir;
  std::string initdatadir;

  const int delta_xy = 120; //There are 120 30 arc-second pieces in one degree
  const double dy    = 6370000.*pi/(180.*delta_xy); //radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.

  double deltat = std::numeric_limits<double>::signaling_NaN();

  void print() const;
};

Parameters LoadParameters(const std::string config_file);

#endif
