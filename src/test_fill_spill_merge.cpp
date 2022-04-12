#include "doctest.h"
#include "irf.cpp"
#include "transient_groundwater.hpp"

#include <fmt/core.h>
#include <richdem/depressions/Zhou2016.hpp>
#include <richdem/terrain_generation.hpp>

#include <random>
#include <sstream>

using namespace richdem;
using namespace richdem::dephier;

#ifdef CODE_COVERAGE
#pragma message \
    "FSM is using a small number of test cases for code coverage estimation. Disable code coverage to enable more extensive tests."
const int number_of_small_tests = 50;
const int number_of_large_tests = 5;
#else
#pragma message \
    "FSM is using a large number of test cases to judge correctness. Enabling code coverage will reduce the number of test cases used."
constexpr int number_of_small_tests = 1000;
constexpr int number_of_large_tests = 200;
#endif

template <class T>
double MaxArrayDiff(const Array2D<T>& a, const Array2D<T>& b) {
  double max_diff = 0;
  for (auto i = a.i0(); i < a.size(); i++) {
    max_diff = std::max(max_diff, (double)std::abs(a(i) - b(i)));
  }
  return max_diff;
}

template <class T>
bool ArrayValuesEqual(const Array2D<T>& a, const Array2D<T>& b) {
  for (auto i = a.i0(); i < a.size(); i++) {
    if (a(i) != b(i))
      return false;
  }
  return true;
}

std::mt19937_64 test_prng;

Array2D<double> random_terrain(std::mt19937_64& gen, const int min_size, const int max_size) {
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  return perlin(size_dist(gen), seed_dist(gen));
}

Array2D<double> random_integer_terrain(std::mt19937_64& gen, const int min_size, const int max_size) {
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  auto dem = perlin(size_dist(gen), seed_dist(gen));
  for (auto i = dem.i0(); i < dem.size(); i++) {
    dem(i) *= 100;
    dem(i) = static_cast<int>(dem(i));
  }

  return dem;
}

TEST_CASE("Determine water level") {
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  const Array2D<double> porosity = {
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
      {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3},
  };

  const std::vector<double> cell_area = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  Array2D<double> wtd = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  };

  const int cx                            = 3;
  const int cy                            = 4;
  const double current_volume             = 5.;
  const double current_area               = 5.;
  const double area_times_elevation_total = 5.;

  SUBCASE("Depression volume exactly equals water volume") {
    const auto water_level = DetermineWaterLevel(
        topo, porosity, cell_area, wtd, 5, cx, cy, current_volume, current_area, area_times_elevation_total);
    CHECK(water_level == 2);  // Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is less than the depression volume") {
    const auto water_level = DetermineWaterLevel(
        topo, porosity, cell_area, wtd, 4, cx, cy, current_volume, current_area, area_times_elevation_total);
    CHECK(water_level == 9 / 5.0);
  }

  SUBCASE("Water volume is greater than the depression volume") {
    const auto water_level = DetermineWaterLevel(
        topo, porosity, cell_area, wtd, 5.5, cx, cy, current_volume, current_area, area_times_elevation_total);
    CHECK(water_level == 2);
  }
}

TEST_CASE("Backfill Depression") {
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9, 6, 6, 6, 6, 6, 6, 6, 6, -9},
      {-9, 6, 1, 6, 1, 1, 6, 1, 6, -9},
      {-9, 6, 1, 6, 1, 3, 6, 1, 6, -9},
      {-9, 6, 1, 6, 2, 1, 4, 1, 6, -9},
      {-9, 6, 1, 6, 1, 1, 6, 1, 6, -9},
      {-9, 6, 1, 6, 6, 6, 6, 1, 6, -9},
      {-9, 6, 1, 1, 1, 1, 1, 1, 6, -9},
      {-9, 6, 6, 6, 6, 6, 6, 6, 6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<double> wtd(topo.width(), topo.height());

  std::vector<flat_c_idx> cells_affected = {
      topo.xyToI(4, 2),
      topo.xyToI(5, 2),
      topo.xyToI(4, 3),
      topo.xyToI(5, 3),
      topo.xyToI(4, 4),
      topo.xyToI(5, 4),
      topo.xyToI(4, 5),
      topo.xyToI(5, 5),
  };

  BackfillDepression(4.0, topo, wtd, cells_affected);

  const Array2D<double> wtd_good = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 3, 3, 0, 0, 0, 0},
      {0, 0, 0, 0, 3, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 2, 3, 0, 0, 0, 0},
      {0, 0, 0, 0, 3, 3, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  };

  CHECK(ArrayValuesEqual(wtd, wtd_good));
}

// TEST_CASE("Add recharge"){
// //double add_recharge(const double deltat, const double my_rech, double my_wtd, const int my_mask, const double
// my_porosity){

//   //positive wtd, positive rech
//   double new_wtd = add_recharge(31536000., 1., 1., 1, 0.4);
//   CHECK(new_wtd==2);

//   //positive wtd, negative rech that still ends positive
//   new_wtd = add_recharge(31536000., -0.5, 1., 1, 0.4);
//   CHECK(new_wtd==0.5);

//   //positive wtd, negative rech that would end negative so gets corrected to 0
//   new_wtd = add_recharge(31536000., -2., 1., 1, 0.4);
//   CHECK(new_wtd==0);

//   //negative wtd, negative rech: wtd should not change
//   new_wtd = add_recharge(31536000., -1., -1., 1, 0.4);
//   CHECK(new_wtd==-1);

//   //negative wtd, positive rech, wtd stays negative
//   new_wtd = add_recharge(31536000., 1., -10., 1, 0.4);
//   CHECK(new_wtd==-7.5);

//   //negative wtd, positive rech, wtd becomes positive
//   new_wtd = add_recharge(31536000., 1., -1., 1, 0.4);
//   CHECK(new_wtd==0.6);
// }

TEST_CASE("Fill a full depression") {
  SubtreeDepressionInfo my_stdi;
  my_stdi.my_labels.emplace(1);
  my_stdi.leaf_label = 1;
  my_stdi.top_label  = 1;

  ArrayPack arp;

  SUBCASE("A very basic depression") {
    arp.topo.resize(10, 10, 2);
    arp.topo.setEdges(-9);

    arp.topo(4, 2) = 1;
    arp.topo(4, 3) = 1;
    arp.topo(4, 4) = 1;
    arp.topo(4, 5) = 1;
    arp.topo(4, 6) = 1;
    arp.topo(3, 6) = 3;

    arp.wtd   = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.label.setEdges(OCEAN);

    arp.final_label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.final_label.setEdges(OCEAN);

    arp.flowdirs  = Array2D<flowdir_t>(arp.topo.width(), arp.topo.height(), NO_FLOW);
    arp.cell_area = std::vector<double>(arp.topo.height(), 1.);

    const Array2D<double> wtd_good = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    auto DH =
        GetDepressionHierarchy<float, Topology::D8>(arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

    FillAFullDepression(my_stdi, DH, arp);

    CHECK(ArrayValuesEqual(arp.wtd, wtd_good));
  }
  SUBCASE("Depression with a little more going on") {
    arp.wtd = Array2D<double>(10, 10, 0);
    arp.topo.resize(arp.wtd.width(), arp.wtd.height(), 4);
    arp.topo.setEdges(-9);

    arp.topo = {
        {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
        {-9, 4, 4, 3, 3, 3, 3, 3, 4, -9},
        {-9, 4, 3, 3, 1, 2, 0, 3, 4, -9},
        {-9, 3, 3, 3, 1, 2, 0, 1, 4, -9},
        {-9, 3, 4, 3, 1, 0, 0, 1, 4, -9},
        {-9, 4, 5, 4, 1, 2, 0, 2, 4, -9},
        {-9, 4, 5, 5, 1, 2, 0, 2, 4, -9},
        {-9, 4, 5, 6, 6, 3, 3, 3, 4, -9},
        {-9, 4, 4, 6, 6, 3, 3, 3, 4, -9},
        {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
    };

    arp.wtd   = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.label.setEdges(OCEAN);

    arp.final_label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.final_label.setEdges(OCEAN);

    arp.flowdirs  = Array2D<flowdir_t>(arp.topo.width(), arp.topo.height(), NO_FLOW);
    arp.cell_area = std::vector<double>(arp.topo.height(), 1.);

    const Array2D<double> wtd_good = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 1, 3, 0, 0, 0},
        {0, 0, 0, 0, 2, 1, 3, 2, 0, 0},
        {0, 0, 0, 0, 2, 3, 3, 2, 0, 0},
        {0, 0, 0, 0, 2, 1, 3, 1, 0, 0},
        {0, 0, 0, 0, 2, 1, 3, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    auto DH =
        GetDepressionHierarchy<float, Topology::D8>(arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

    FillAFullDepression(my_stdi, DH, arp);

    CHECK(ArrayValuesEqual(arp.wtd, wtd_good));
  }
}

void CheckMassLoss(const int count, const int min_size, const int max_size) {
  for (int test_num = 0; test_num < count; test_num++) {
    std::stringstream oss;

    ArrayPack arp;
    Parameters params;

    Array2D<double> dem;

    {
      oss << test_prng;
      dem = random_terrain(test_prng, min_size, max_size);
      std::cerr << "checking mass loss #" << test_num << std::endl;
    }

    arp.topo  = dem;
    arp.label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.label.setEdges(OCEAN);
    arp.final_label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.final_label.setEdges(OCEAN);
    arp.flowdirs        = Array2D<flowdir_t>(arp.topo.width(), arp.topo.height(), NO_FLOW);
    arp.wtd             = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.wtd_T           = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.wtd_T_iteration = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.original_wtd    = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);

    arp.precip                = Array2D<double>(arp.topo.width(), arp.topo.height(), 1.5);
    arp.starting_evap         = Array2D<double>(arp.topo.width(), arp.topo.height(), 1.0);
    arp.open_water_evap       = Array2D<double>(arp.topo.width(), arp.topo.height(), 1.3);
    arp.ksat                  = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.00001);
    arp.fdepth                = Array2D<double>(arp.topo.width(), arp.topo.height(), 50.);
    arp.porosity              = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.3);
    arp.effective_storativity = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.3);
    arp.rech                  = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.5);
    arp.transmissivity        = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);

    arp.runoff         = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.scalar_array_y = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);
    arp.scalar_array_x = Array2D<double>(arp.topo.width(), arp.topo.height(), 0.);

    arp.cell_area.resize(arp.topo.height());
    arp.cellsize_e_w_metres.resize(arp.topo.height());
    arp.land_mask = Array2D<uint8_t>(arp.topo.width(), arp.topo.height(), 1);

    for (unsigned int j = 0; j < arp.cell_area.size(); j++) {
      arp.cell_area[j]           = 100;
      arp.cellsize_e_w_metres[j] = 10;
    }

    params.ncells_x            = arp.topo.width();
    params.ncells_y            = arp.topo.height();
    params.cellsize_n_s_metres = 10;
    params.infiltration_on     = false;

    for (int y = 0; y < params.ncells_y; y++) {
      for (int x = 0; x < params.ncells_x; x++) {
        arp.topo(x, y)                  = dem(x, y);
        arp.precip(x, y)                = 1.5;
        arp.starting_evap(x, y)         = 1.;
        arp.open_water_evap(x, y)       = 1.3;
        arp.ksat(x, y)                  = 0.00001;
        arp.fdepth(x, y)                = 50.;
        arp.porosity(x, y)              = 0.3;
        arp.effective_storativity(x, y) = 0.3;
        arp.rech(x, y)                  = 0.5;
        arp.land_mask(x, y)             = 1;
        if (arp.topo(x, y) == 0) {
          arp.label(x, y)       = OCEAN;
          arp.final_label(x, y) = OCEAN;
          arp.land_mask(x, y)   = 0;
        }
      }
    }

    arp.topo.setEdges(0);
    arp.land_mask.setEdges(0);

    params.deltat = 10000;

    auto deps =
        GetDepressionHierarchy<float, Topology::D8>(arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);
    std::cerr << "DH done" << std::endl;

    double change_in_rech        = 0.;
    double change_in_ocean_loss  = 0;
    double change_in_water_table = 0;
    double rech_old              = 0;
    double ocean_loss_old        = 0;
    double water_table_old       = 0;
    double mass_loss             = 0.;

    for (int reps = 0; reps < 10; reps++) {
      FanDarcyGroundwater::update(params, arp);
      std::cout << "GW done" << std::endl;
      dh::FillSpillMerge(params, deps, arp);
      std::cout << "FSM done" << std::endl;

      if (params.evap_mode) {
        std::cout << "updating the evaporation field" << std::endl;
#pragma omp parallel for default(none) shared(arp)
        for (unsigned int i = 0; i < arp.topo.size(); i++) {
          if (arp.wtd(i) > 0) {  // if there is surface water present
            arp.rech(i) = arp.precip(i) - arp.open_water_evap(i);
          } else {  // water table is below the surface
            arp.rech(i) = arp.precip(i) - arp.starting_evap(i);
            if (arp.rech(i) < 0) {  // Recharge is always positive.
              arp.rech(i) = 0.;
            }
          }
        }
      }

      double wtd_sum = 0.0;
      for (int y = 0; y < params.ncells_y; y++) {
        for (int x = 0; x < params.ncells_x; x++) {
          if (arp.wtd(x, y) > 0)
            wtd_sum += arp.wtd(x, y) * arp.cell_area[y];
          else
            wtd_sum += arp.wtd(x, y) * arp.porosity(x, y) * arp.cell_area[y];
        }
      }

      change_in_rech        = params.total_added_recharge - rech_old;
      change_in_ocean_loss  = params.total_loss_to_ocean - ocean_loss_old;
      change_in_water_table = wtd_sum - water_table_old;

      mass_loss = (change_in_rech - change_in_water_table - change_in_ocean_loss) / change_in_rech * 100.;

      std::cout << "change in rech is " << change_in_rech << " total_added_recharge " << params.total_added_recharge
                << " rech old " << rech_old << std::endl;
      std::cout << "change in ocean loss " << change_in_ocean_loss << " total loss to ocean "
                << params.total_loss_to_ocean << " ocean loss old " << ocean_loss_old << std::endl;
      std::cout << "change in water table " << change_in_water_table << " wtd sum " << wtd_sum << " water table old "
                << water_table_old << std::endl;
      std::cout << "mass loss as a percentage of recharge is " << mass_loss << std::endl;

      water_table_old = wtd_sum;
      ocean_loss_old  = params.total_loss_to_ocean;
      rech_old        = params.total_added_recharge;
    }

    CHECK_MESSAGE(mass_loss < 100, "failed mass loss check");
  }
}

TEST_CASE("Check mass loss for random cases") {
  CheckMassLoss(number_of_small_tests, 10, 30);
  CheckMassLoss(number_of_large_tests, 100, 300);
}

// This is hard to make work correctly: we need to know which depressions are contained within which for the function.
// void RandomizedFullDepressionVsPriorityFlood(const int count, const int min_size, const int max_size){
//  for(int i=0;i<count;i++){
//  std::stringstream oss;
//
//  ArrayPack arp;
//
//    Array2D<double> topo;
//
//  {
//    oss<<gen;
//    topo = random_terrain(gen, min_size, max_size);
//    std::cerr<<"Randomized Full Depression vs Priority-Flood #"<<i<<std::endl;
//  }
////topo.printAll("topo before");
//
// arp.topo.resize(topo.width(),topo.height(),0);
//    for (int i=0;i<arp.topo.size();++i)
//        arp.topo(i)=static_cast<float> (topo(i));
////arp.topo.printAll("topo arp");
//
//
//    arp.topo.setEdges(-9);
//    arp.label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
//    arp.label.setEdges(OCEAN);
//    arp.final_label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
//    arp.final_label.setEdges(OCEAN);
//    arp.flowdirs = Array2D<flowdir_t>  (arp.topo.width(), arp.topo.height(), NO_FLOW);
//    arp.cell_area = std::vector<double> (arp.topo.height(), 1.);
//    arp.wtd = Array2D<double> (arp.topo.width(), arp.topo.height(), 0. );
// std::cerr<<"ive set all of those"<<std::endl;
//
//    auto DH = GetDepressionHierarchy<float,Topology::D8>(arp.topo, arp.cell_area, arp.label, arp.final_label,
//    arp.flowdirs); std::cerr<<"DH done"<<std::endl;
//
//
//    for(i=1;i<arp.label.size();i++){
//      SubtreeDepressionInfo my_stdi;
//      my_stdi.my_labels.emplace(arp.label(i));
//      my_stdi.top_label = arp.label(i);
//      FillAFullDepression(my_stdi, DH, arp);
// std::cerr<<std::setprecision(15)<<"and set the labels, "<<arp.label(i)<<" my out cell is
// "<<DH.at(my_stdi.top_label).out_elev<<std::endl;
//    	}
//
//
// arp.label.printAll("label");
// arp.topo.printAll("topo");
//
//
//
// std::cerr<<"full depression done"<<std::endl;
// std::cerr<<std::setprecision(15)<<"test mine before adding "<<arp.topo(5,5)<<std::endl;
//
//    //show the filled topo
//    for(auto i=arp.topo.i0(); i<arp.topo.size(); i++){
//      if(!arp.topo.isNoData(i))
//        arp.topo(i) += arp.wtd(i);
//    }
//
//
//
//
//
//    auto comparison_dem = arp.topo;
//    PriorityFlood_Zhou2016(comparison_dem);
//
// arp.topo.printAll("wtd");
// comparison_dem.printAll("priority flood dem");
//
// std::cerr<<std::setprecision(15)<<"test mine "<<arp.topo(5,5)<<std::endl;
// std::cerr<<std::setprecision(15)<<"test priority "<<comparison_dem(5,5)<<std::endl;
// std::cerr<<std::setprecision(15)<<"test wtd "<<arp.wtd(5,5)<<std::endl;
//
//
//    CHECK_MESSAGE(MaxArrayDiff(comparison_dem,arp.topo)<1e-6, "Randomized full depression vs Priority-Flood failed
//    with width = "+std::to_string(arp.topo.width())+" height =
//    "+std::to_string(arp.topo.height())+" state = " + oss.str());
//  }
//}
//
//
// TEST_CASE("Randomized Full depression vs Priority-Flood"){
//  RandomizedFullDepressionVsPriorityFlood(number_of_small_tests,  10,  30);
//  RandomizedFullDepressionVsPriorityFlood(number_of_large_tests, 100, 300);
//}

void RandomizedHeavyFloodingVsPriorityFlood(const int count, const int min_size, const int max_size) {
#pragma omp parallel for
  for (int test_num = 0; test_num < count; test_num++) {
    std::stringstream oss;

    ArrayPack arp;
    Parameters params;
    params.textfilename = "test.txt";

    Array2D<double> dem;

#pragma omp critical
    {
      oss << test_prng;
      dem = random_terrain(test_prng, min_size, max_size);
      std::cerr << "Randomized Heavy Flooding vs Priority-Flood #" << test_num << std::endl;
    }

    arp.topo        = dem;
    arp.label       = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.final_label = Array2D<dh_label_t>(arp.topo.width(), arp.topo.height(), NO_DEP);
    arp.flowdirs    = Array2D<flowdir_t>(arp.topo.width(), arp.topo.height(), NO_FLOW);

    arp.cell_area = std::vector<double>(arp.topo.height(), 1.);

    arp.topo.setEdges(-1);
    arp.label.setEdges(OCEAN);

    auto deps =
        GetDepressionHierarchy<double, Topology::D8>(arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

    // wtd with a *lot* of initial surface water
    arp.wtd.resize(arp.topo.width(), arp.topo.height());
    arp.wtd.setAll(0);

    arp.runoff.resize(arp.topo.width(), arp.topo.height());
    arp.runoff.setAll(0);

    arp.porosity.resize(arp.topo.width(), arp.topo.height());
    arp.porosity.setAll(1);

    try {
      FillSpillMerge(params, deps, arp);
    } catch (const std::exception& e) {
      std::cerr << fmt::format(
                       "FillSpillMerge failed because of \"{}\" with width = {}, height = {}, state = {}",
                       e.what(),
                       arp.topo.width(),
                       arp.topo.height(),
                       oss.str())
                << std::endl;
      throw e;
    }

    for (auto i = arp.topo.i0(); i < arp.topo.size(); i++) {
      if (!arp.topo.isNoData(i))
        arp.topo(i) += arp.wtd(i);
    }

    auto comparison_dem = arp.topo;
    PriorityFlood_Zhou2016(comparison_dem);

    CHECK_MESSAGE(
        MaxArrayDiff(comparison_dem, arp.topo) < 1e-6,
        fmt::format(
            "Randomized Heavy Flooding vs Priority-Flood failed with width = {}, height = {}, state = {}",
            arp.topo.width(),
            arp.topo.height(),
            oss.str()));
  }
}

TEST_CASE("Randomized Heavy Flooding vs Priority-Flood") {
  RandomizedHeavyFloodingVsPriorityFlood(number_of_small_tests, 10, 30);
  RandomizedHeavyFloodingVsPriorityFlood(number_of_large_tests, 100, 300);
}
