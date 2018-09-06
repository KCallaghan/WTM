#include "Array2D.hpp"
#include "dephier.hpp"
#include "DisjointDenseIntSet.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  const std::string in_name   = argv[1];
  const std::string out_name  = argv[2];
  const std::string out_graph = argv[3];

  Array2D<float> dem(in_name,"value");   //Recharge (Percipitation minus Evapotranspiration)

  //Initialize labels to indicate that none of the cells are part of depressions
  Array2D<label_t> label   (dem.width(),dem.height(),NO_DEP);

  //Initialize flow directions to indicate that none of the cells flow anywhere
  Array2D<flowdir_t> flowdirs(dem.width(),dem.height(),NO_FLOW);

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(int i=0;i<label.size();i++)
    if(dem(i)==0)
      label(i) = OCEAN;

  //Label all the depressions and get the hierarchy connecting them.
  const auto deps = GetDepressionHierarchy<float,Topology::D8>(dem, label, flowdirs);

  //TODO: Remove. For viewing test cases.
  if(label.width()<1000){
    for(int y=0;y<label.height();y++){
      for(int x=0;x<label.width();x++)
        std::cout<<std::setw(3)<<label(x,y)<<" ";
      std::cout<<std::endl;
    }

    //GraphViz dot-style output for drawing depression hierarchy graphs.
    std::ofstream fgraph(out_graph);
    fgraph<<"digraph {\n";
    for(unsigned int i=0;i<deps.size();i++)
      fgraph<<i<<" -> "<<deps[i].parent<<";\n";
    fgraph<<"}\n";
  }

  SaveAsNetCDF(dem,out_name+"-dem.nc","value");
  SaveAsNetCDF(label,out_name+"-labels_raw.nc","value");

  LastLayer(label, dem, deps);

  SaveAsNetCDF(label,out_name+"-labels_proc.nc","value");

  return 0;
}
