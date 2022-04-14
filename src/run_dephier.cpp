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

  constexpr float earth_radius = 6371000.;  // metres
  // Radius * Pi = Distance from N to S pole
  // Distance / 180 = Meters / degree latitude
  const auto meters_per_degree = earth_radius * DEG_TO_RAD;

  // distance between lines of latitude is a constant.
  params.cellsize_n_s_metres = meters_per_degree / params.cells_per_degree;

  // initialise some arrays
  // size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cellsize_e_w_metres.resize(params.ncells_y);
  // cell area (metres squared)
  arp.cell_area.resize(params.ncells_y);

  for (int32_t j = 0; j < params.ncells_y; j++) {
    // southern edge of the domain in degrees, plus the number of cells up
    // from this location/the number of cells per degree, converted to radians.
    // latitude at the southern edge of a cell (subtract half a cell):
    const double latitude_radians_S = ((float(j)) / params.cells_per_degree + params.southern_edge) * DEG_TO_RAD;
    // latitude at the northern edge of a cell (add half a cell):
    const double latitude_radians_N = ((float(j) + 1.) / params.cells_per_degree + params.southern_edge) * DEG_TO_RAD;

    // distance at the northern edge of the cell for the given latitude:
    double cellsize_e_w_metres_N = params.cellsize_n_s_metres * std::cos(latitude_radians_N);
    // distance at the southern edge of the cell for the given latitude:
    double cellsize_e_w_metres_S = params.cellsize_n_s_metres * std::cos(latitude_radians_S);

    // distance between lines of longitude varies with latitude.
    // This is the distance at the centre of a cell for a given latitude:
    arp.cellsize_e_w_metres[j] = (cellsize_e_w_metres_N + cellsize_e_w_metres_S) / 2.;

    // cell area computed as a trapezoid, using unchanging north-south distance,
    // and east-west distances at the northern and southern edges of the cell:
    arp.cell_area[j] = params.cellsize_n_s_metres * arp.cellsize_e_w_metres[j];
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
