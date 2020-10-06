#include "doctest.h"
#include "irf.cpp"

#include <richdem/terrain_generation.hpp>
#include <random>
#include <sstream>

using namespace richdem;
using namespace richdem::dephier;


#ifdef CODE_COVERAGE
  #pragma message "FSM is using a small number of test cases for code coverage estimation. Disable code coverage to enable more extensive tests."
  const int number_of_small_tests = 50;
  const int number_of_large_tests = 5;
#else
  #pragma message "FSM is using a large number of test cases to judge correctness. Enabling code coverage will reduce the number of test cases used."
  const int number_of_small_tests = 6000;
  const int number_of_large_tests = 500;
#endif



template<class T>
double MaxArrayDiff(const Array2D<T> &a, const Array2D<T> &b){
  double max_diff = 0;
  for(auto i=a.i0();i<a.size();i++){
    max_diff = std::max(max_diff, (double)std::abs(a(i)-b(i)));
  }
  return max_diff;
}

template<class T>
bool ArrayValuesEqual(const Array2D<T> &a, const Array2D<T> &b){
  for(auto i=a.i0();i<a.size();i++){
    if(a(i)!=b(i))
      return false;
  }
  return true;
}



std::mt19937_64 gen;

Array2D<double> random_terrain(std::mt19937_64 &gen, const int min_size, const int max_size){
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  return perlin(size_dist(gen), seed_dist(gen));
}

Array2D<double> random_integer_terrain(std::mt19937_64 &gen, const int min_size, const int max_size){
  static std::uniform_int_distribution<uint32_t> seed_dist;

  std::uniform_int_distribution<int> size_dist(min_size, max_size);

  auto dem = perlin(size_dist(gen), seed_dist(gen));
  for(auto i=dem.i0();i<dem.size();i++){
    dem(i) *= 100;
    dem(i) = static_cast<int>(dem(i));
  }

  return dem;
}



TEST_CASE("Determine water level"){

  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  2,  2,  2,  2,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  1,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  1,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  1,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  1,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  1,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  2,  2,  2,  2,  2, -9},
      {-9,  2,  2,  2,  2,  2,  2,  2,  2, -9},
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

  const std::vector<double> cell_area = {1,1,1,1,1,1,1,1,1,1};


  Array2D<double> wtd = {
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, -1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0,  0, 0, 0, 0, 0, 0, 0},
  };


  const int cx = 3;
  const int cy = 4;
  const double current_volume = 5.;
  const double current_area = 5.;
  const double area_times_elevation_total = 5.;


  SUBCASE("Depression volume exactly equals water volume"){
    const auto water_level = DetermineWaterLevel(topo,porosity,cell_area,wtd,5,cx,cy,current_volume,current_area,area_times_elevation_total);
    CHECK(water_level==2); //Water elevation equals the sill elevation
  }

  SUBCASE("Water volume is less than the depression volume"){
    const auto water_level = DetermineWaterLevel(topo,porosity,cell_area,wtd,4,cx,cy,current_volume,current_area,area_times_elevation_total);
    CHECK(water_level==9/5.0);
  }

  SUBCASE("Water volume is greater than the depression volume"){
    const auto water_level = DetermineWaterLevel(topo,porosity,cell_area,wtd,5.5,cx,cy,current_volume,current_area,area_times_elevation_total);
    CHECK(water_level==2);
  }
}




TEST_CASE("Backfill Depression"){
  const Array2D<double> topo = {
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  1,  3,  6,  1,  6, -9},
      {-9,  6,  1,  6,  2,  1,  4,  1,  6, -9},
      {-9,  6,  1,  6,  1,  1,  6,  1,  6, -9},
      {-9,  6,  1,  6,  6,  6,  6,  1,  6, -9},
      {-9,  6,  1,  1,  1,  1,  1,  1,  6, -9},
      {-9,  6,  6,  6,  6,  6,  6,  6,  6, -9},
      {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9},
  };

  Array2D<double> wtd(topo.width(), topo.height());

  std::vector<flat_c_idx> cells_affected = {
    topo.xyToI(4,2), topo.xyToI(5,2),
    topo.xyToI(4,3), topo.xyToI(5,3),
    topo.xyToI(4,4), topo.xyToI(5,4),
    topo.xyToI(4,5), topo.xyToI(5,5),
  };

  BackfillDepression(4.0, topo, wtd, cells_affected);

  const Array2D<double> wtd_good = {
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  1,  0,  0,  0,  0},
      { 0,  0,  0,  0,  2,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  3,  3,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
  };

  CHECK(ArrayValuesEqual(wtd,wtd_good));
}











void RandomizedHeavyFloodingVsPriorityFlood(const int count, const int min_size, const int max_size){
  #pragma omp parallel for
  for(int i=0;i<count;i++){
    std::stringstream oss;

    ArrayPack arp;
    Parameters params;

    Array2D<double> dem;

    #pragma omp critical
    {
      oss<<gen;
      dem = random_terrain(gen, min_size, max_size);
      std::cerr<<"Randomized Heavy Flooding vs Priority-Flood #"<<i<<std::endl;
    }

    arp.topo = dem;
    arp.label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
    arp.final_label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
    arp.flowdirs = Array2D<flowdir_t>  (arp.topo.width(), arp.topo.height(), NO_FLOW);

    arp.cell_area = std::vector<double> (arp.topo.height(), 1.);


    arp.topo.setEdges(-1);
    arp.label.setEdges(OCEAN);

    auto deps = GetDepressionHierarchy<double,Topology::D8>(arp.topo, arp.cell_area,arp.label,arp.final_label, arp.flowdirs);

    //wtd with a *lot* of initial surface water
    Array2D<double> wtd(arp.topo.width(), arp.topo.height(), 100);

    try {
      FillSpillMerge(params, deps, arp);
    } catch (const std::exception &e) {
      std::cerr<<"FillSpillMerge failed because of \""<<e.what()<<"\" with width = "<<arp.topo.width()<<" height = "<<arp.topo.height()<<" state = "<<oss.str()<<std::endl;
      throw e;
    }

    for(auto i=arp.topo.i0(); i<arp.topo.size(); i++){
      if(!arp.topo.isNoData(i))
        arp.topo(i) += wtd(i);
    }

    auto comparison_dem = arp.topo;
    PriorityFlood_Zhou2016(comparison_dem);

    CHECK_MESSAGE(MaxArrayDiff(comparison_dem,arp.topo)<1e-6, "Randomized Heavy Flooding vs Priority-Flood failed with width = "+std::to_string(arp.topo.width())+" height = "+std::to_string(arp.topo.height())+" state = " + oss.str());
  }
}

TEST_CASE("Randomized Heavy Flooding vs Priority-Flood"){
  RandomizedHeavyFloodingVsPriorityFlood(number_of_small_tests,  10,  30);
  RandomizedHeavyFloodingVsPriorityFlood(number_of_large_tests, 100, 300);
}














