#include "ArrayPack.hpp"
#include "fill_spill_merge.hpp"
#include "parameters.hpp"
#include "transient_groundwater.hpp"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/timer.hpp>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace rd = richdem;
namespace dh = richdem::dephier;

constexpr double DEG_TO_RAD = M_PI / 180.0;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Syntax: " << argv[0] << " <Configuration File>" << std::endl;
    return -1;
  }

  ArrayPack arp;
  std::cerr << "Argv" << argv << std::endl;
  Parameters params(argv[1]);

  // load in the data files: topography and mask.

  std::cout << params.surfdatadir << std::endl;

  arp.topo = rd::Array2D<float>(params.surfdatadir + params.region + params.time_start + "_topography.tif");

  arp.land_mask = rd::Array2D<float>(params.surfdatadir + params.region + params.time_start + "_mask.tif");

  //  arp.topo          = LoadData<float>(params.surfdatadir + params.region
  //  + params.time_start + "_topo.nc",    "value");
  //  arp.land_mask     = LoadData<uint8_t>(params.surfdatadir + params.region +
  //  params.time_start + "_mask.nc",   "value");

  // width and height in number of cells in the array
  params.ncells_x = arp.topo.width();
  params.ncells_y = arp.topo.height();

  // initialise the label and flow direction arrays:
  rd::Array2D<dh::dh_label_t> label(params.ncells_x, params.ncells_y, dh::NO_DEP);  // No cells are part of a depression

  rd::Array2D<dh::dh_label_t> final_label(
      params.ncells_x, params.ncells_y, dh::NO_DEP);  // No cells are part of a depression

  rd::Array2D<rd::flowdir_t> flowdirs(params.ncells_x, params.ncells_y, rd::NO_FLOW);  // No cells flow anywhere

  const float earth_radius = 6371000.;  // metres

  // distance between lines of latitude is a constant.
  params.cellsize_n_s_metres = (earth_radius * DEG_TO_RAD) / params.cells_per_degree;

  // initialise some arrays
  arp.latitude_radians.resize(params.ncells_y);
  // the latitude of each row of cells
  arp.cellsize_e_w_metres.resize(params.ncells_y);
  // size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cell_area.resize(params.ncells_y);
  // cell area (metres squared)

  for (unsigned int j = 0; j < arp.latitude_radians.size(); j++) {
    // latitude at the centre of a cell:
    arp.latitude_radians[j] = (float(j) / params.cells_per_degree + params.southern_edge) * DEG_TO_RAD;
    // southern edge of the domain in degrees, plus the number of cells up
    // from this location/the number of cells per degree, converted to radians.

    // cells_per_degree = 120, there are this many 30 arc-second pieces in
    // one degree. (or however many pieces per degree the user has)
    // j/cells_per_degree gives the number  of degrees up from the southern edge,
    // add southern_edge since the southern edge may not be at 0 latitude.
    //*pi/180 to convert to radians.
    // latitude_radians is now the latitude in radians.

    // latitude at the southern edge of a cell (subtract half a cell):
    double latitude_radians_S = ((float(j) - 0.5) / params.cells_per_degree + params.southern_edge) * DEG_TO_RAD;
    // latitude at the northern edge of a cell (add half a cell):
    double latitude_radians_N = ((float(j) + 0.5) / params.cells_per_degree + params.southern_edge) * DEG_TO_RAD;

    // distance between lines of longitude varies with latitude.
    // This is the distance at the centre of a cell for a given latitude:
    arp.cellsize_e_w_metres[j] =
        earth_radius * std::cos(arp.latitude_radians[j]) * DEG_TO_RAD / params.cells_per_degree;

    // distance at the northern edge of the cell for the given latitude:
    double cellsize_e_w_metres_N = earth_radius * std::cos(latitude_radians_N) * DEG_TO_RAD / params.cells_per_degree;
    // distance at the southern edge of the cell for the given latitude:
    double cellsize_e_w_metres_S = earth_radius * std::cos(latitude_radians_S) * DEG_TO_RAD / params.cells_per_degree;

    // cell area computed as a trapezoid, using unchanging north-south distance,
    // and east-west distances at the northern and southern edges of the cell:
    arp.cell_area[j] = params.cellsize_n_s_metres * (cellsize_e_w_metres_N + cellsize_e_w_metres_S) / 2;

    // TODO(kcallaghan): see which, if any, arrays can be cleared after use to free up memory.
    // How to do this?
  }

  // Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
#pragma omp parallel for default(none) shared(arp, final_label, label)
  for (unsigned int i = 0; i < label.size(); i++) {
    if (arp.land_mask(i) == 0) {
      label(i)       = dh::OCEAN;
      final_label(i) = dh::OCEAN;
    }
  }

  // Generate flow directions, label all the depressions, and get the hierarchy
  // connecting them
  std::cout << "going to do depression hierarchy" << std::endl;
  auto deps =
      dh::GetDepressionHierarchy<float, rd::Topology::D8>(arp.topo, arp.cell_area, label, final_label, flowdirs);

  // We are finished, save the result.
  std::cout << "done with processing" << std::endl;

  //  label.saveGDAL("label.tif");
  //  final_label.saveGDAL("final_label.tif");

  //  SaveAsNetCDF(label,"label.nc","value");
  //  SaveAsNetCDF(final_label,"final_label.nc","value");

  return 0;
}
