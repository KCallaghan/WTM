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



TEST_CASE("Add recharge"){
//double add_recharge(const double deltat, const double my_rech, double my_wtd, const int my_mask, const double my_porosity){

  //positive wtd, positive rech
  double new_wtd = add_recharge(31536000., 1., 1., 1, 0.4);
  CHECK(new_wtd==2);

  //positive wtd, negative rech that still ends positive
  new_wtd = add_recharge(31536000., -0.5, 1., 1, 0.4);
  CHECK(new_wtd==0.5);

  //positive wtd, negative rech that would end negative so gets corrected to 0
  new_wtd = add_recharge(31536000., -2., 1., 1, 0.4);
  CHECK(new_wtd==0);

  //negative wtd, negative rech: wtd should not change
  new_wtd = add_recharge(31536000., -1., -1., 1, 0.4);
  CHECK(new_wtd==-1);

  //negative wtd, positive rech, wtd stays negative
  new_wtd = add_recharge(31536000., 1., -10., 1, 0.4);
  CHECK(new_wtd==-7.5);

  //negative wtd, positive rech, wtd becomes positive
  new_wtd = add_recharge(31536000., 1., -1., 1, 0.4);
  CHECK(new_wtd==0.6);

}


TEST_CASE("Fill a full depression"){
  SubtreeDepressionInfo my_stdi;
  my_stdi.my_labels.emplace(1);
  my_stdi.leaf_label = 1;
  my_stdi.top_label = 1;

  ArrayPack arp;

  arp.wtd = Array2D<double>(10,10,0);
  arp.topo.resize(arp.wtd.width(),arp.wtd.height(),2);
  arp.topo.setEdges(-9);
  arp.topo(4,2) = 1;
  arp.topo(4,3) = 1;
  arp.topo(4,4) = 1;
  arp.topo(4,5) = 1;
  arp.topo(4,6) = 1;



  std::cout<<"topo"<<std::endl;
  arp.topo.printAll();

  const Array2D<double> dem = {
      {-9,-9,-9,-9,-9,-9,-9,-9,-9, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 1, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9, 2, 2, 2, 2, 2, 2, 2, 2, -9},
      {-9,-9,-9,-9,-9,-9,-9,-9,-9, -9},
  };

 //not sure why passing dem to GetDepressionHierarchy below works, and passing arp.topo doesn't?? So weird
  arp.wtd = Array2D<double> (arp.topo.width(), arp.topo.height(), 0. );
  arp.label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
  arp.label.setEdges(OCEAN);

  arp.final_label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
  arp.final_label.setEdges(OCEAN);

  arp.flowdirs = Array2D<flowdir_t>  (arp.topo.width(), arp.topo.height(), NO_FLOW);
  arp.cell_area = std::vector<double> (arp.topo.height(), 1.);

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

  auto DH = GetDepressionHierarchy<double,Topology::D8>(dem, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

  FillAFullDepression(my_stdi, DH, arp);

std::cout<<"labels"<<std::endl;
  arp.label.printAll();

std::cout<<"wtd"<<std::endl;
  arp.wtd.printAll();




  CHECK(ArrayValuesEqual(arp.wtd,wtd_good));

}





//void RandomizedHeavyFloodingVsPriorityFlood(const int count, const int min_size, const int max_size){
//  #pragma omp parallel for
//  for(int i=0;i<count;i++){
//    std::stringstream oss;
//
//    ArrayPack arp;
//    Parameters params;
//    params.textfilename = "test.txt";
//
//    Array2D<double> dem;
//
//    #pragma omp critical
//    {
//      oss<<gen;
//      dem = random_terrain(gen, min_size, max_size);
//      std::cerr<<"Randomized Heavy Flooding vs Priority-Flood #"<<i<<std::endl;
//    }
//
//    arp.topo = dem;
//    arp.label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
//    arp.final_label = Array2D<dh_label_t> (arp.topo.width(), arp.topo.height(), NO_DEP );
//    arp.flowdirs = Array2D<flowdir_t>  (arp.topo.width(), arp.topo.height(), NO_FLOW);
//
//    arp.cell_area = std::vector<double> (arp.topo.height(), 1.);
//
//
//    arp.topo.setEdges(-1);
//    arp.label.setEdges(OCEAN);
//
//    auto deps = GetDepressionHierarchy<double,Topology::D8>(arp.topo, arp.cell_area,arp.label,arp.final_label, arp.flowdirs);
//
//    //wtd with a *lot* of initial surface water
//    Array2D<double> wtd(arp.topo.width(), arp.topo.height(), 100);
//
//    try {
//      FillSpillMerge(params, deps, arp);
//    } catch (const std::exception &e) {
//      std::cerr<<"FillSpillMerge failed because of \""<<e.what()<<"\" with width = "<<arp.topo.width()<<" height = "<<arp.topo.height()<<" state = "<<oss.str()<<std::endl;
//      throw e;
//    }
//
//    for(auto i=arp.topo.i0(); i<arp.topo.size(); i++){
//      if(!arp.topo.isNoData(i))
//        arp.topo(i) += wtd(i);
//    }
//
//    auto comparison_dem = arp.topo;
//    PriorityFlood_Zhou2016(comparison_dem);
//
//    CHECK_MESSAGE(MaxArrayDiff(comparison_dem,arp.topo)<1e-6, "Randomized Heavy Flooding vs Priority-Flood failed with width = "+std::to_string(arp.topo.width())+" height = "+std::to_string(arp.topo.height())+" state = " + oss.str());
//  }
//}
//
//TEST_CASE("Randomized Heavy Flooding vs Priority-Flood"){
//  RandomizedHeavyFloodingVsPriorityFlood(number_of_small_tests,  10,  30);
//  RandomizedHeavyFloodingVsPriorityFlood(number_of_large_tests, 100, 300);
//}














