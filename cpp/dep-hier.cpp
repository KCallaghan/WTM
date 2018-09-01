//#include "Array2D.hpp"
#include <richdem/common/Array2D.hpp>
#include <richdem/common/grid_cell.hpp>
#include "DisjointDenseIntSet.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

constexpr double SQ2 = std::sqrt(2.0);

using namespace richdem;

//1 2 3
//0   4
//7 6 5
//                      0  1  2  3 4 5 6  7
const int dx8[8]       = {-1,-1, 0, 1,1,1,0,-1};
const int dy8[8]       = {0, -1,-1,-1,0,1,1, 1};
const double dr8[8]    = {1,SQ2,1,SQ2,1,SQ2,1,SQ2};
const int d8inverse[8] = {4,  5, 6, 7,0,1,2, 3};

//  1
//0   2
//  3
//                      0  1 2 3
const int dx4[8]       = {-1, 0,1,0};
const int dy4[8]       = { 0,-1,0,1};
const double dr4[4]    = { 1, 1,1,1};
const int d4inverse[4] = { 2, 3,0,1};

const int    *const mdx      = dx8;
const int    *const mdy      = dy8;
const int    *const dinverse = d8inverse;
const double *const mdr      = dr8;
const int neighbours         = 8;


const float  OCEAN_LEVEL = 0;
const int8_t NO_FLOW     = -1;

class mGridCellZk {
 public:
  int x, y;
  float z;
  mGridCellZk(const int x0, const int y0, const float z0){
    x = x0;
    y = y0;
    z = z0;
  }
  bool operator>(const mGridCellZk& a) const {
    return z>a.z; //Less than sorts the queue in reverse
  }
};

//Identifier for cells which are not part of a depression
const int NO_DEP = -1;
const int OCEAN  = 0;

typedef int label_t;

class Depression {
 public:
  label_t pit_cell = -1;
  label_t out_cell = -1;
  label_t parent   = -1;
  double  pit_elev = std::numeric_limits<double>::infinity();
  double  out_elev = std::numeric_limits<double>::infinity();  
};

class Outlet {
 public:
  label_t depa;
  label_t depb;
  label_t out_cell = -1;
  double  out_elev = std::numeric_limits<double>::infinity();
  Outlet(label_t depa0, label_t depb0, label_t out_cell0, double out_elev0){
    depa     = depa0;
    depb     = depb0;
    out_cell = out_cell0;
    out_elev = out_elev0;
  }
  bool operator==(const Outlet &o) const {
    return (depa==o.depa && depb==o.depb) || (depa==o.depb && depb==o.depa);
  }
};

struct OutletHash {
  std::size_t operator()(const Outlet &out) const {
    return out.depa ^ out.depb;
  }
};

int ModFloor(int a, int n) {
  return ((a % n) + n) % n;
}

template<class T>
void WrapCoordinates(const Array2D<T> &arr, int &nx, int &ny){
  nx = ModFloor(nx,arr.width());
}



std::vector<Depression> GetDepressionHierarchy(const Array2D<float> &topo, Array2D<int> &label){
  //Depressions are identified by a number [0,*). This vector holds the
  //depressions.
  std::vector<Depression> depressions;

  std::unordered_set<Outlet, OutletHash> outlet_database;

  //The priority queue ensures that cells are visited in order from lowest to
  //highest. Cells of equal elevation that are added to the queue later are
  //processed before cells that were added earlier.
  GridCellZk_high_pq<float> pq;

  //Ensure that all cells are initially labeled as being not part of a
  //depression.
  label.setAll(NO_DEP);

  //We start by adding all the edge cells to the priority queue. We label these
  //edge cells as ocean. Note that if a landmass were to intersect the edges of
  //the DEM, this may lead to incorrect results. The user should find a way of
  //identifying what they consider to be the ocean/base-level (if any) and
  //modify this portion of the code to identify it.
  // for(int y=0;y<topo.height();y++){
  //   pq.emplace(0,               y,topo(0,             y)); 
  //   pq.emplace(topo.width()-1,  y,topo(topo.width()-1,y)); 
  //   label(0,               y) = OCEAN;
  //   label(topo.width()-1,  y) = OCEAN;    
  // }
  // for(int x=0;x<topo.width();x++){
  //   pq.emplace(x, 0,               topo(x,               0)); 
  //   pq.emplace(x, topo.height()-1, topo(x, topo.height()-1)); 
  //   label(x, 0              ) = OCEAN;
  //   label(x, topo.height()-1) = OCEAN;    
  // }
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(topo(x,y)!=0)
      continue;
   pq.emplace(x,y,topo(x,y));
   label(x, y) = OCEAN;    
  }

  //Build the ocean depression
  {
    auto &oceandep    = depressions.emplace_back();
    oceandep.pit_cell = -1;
    oceandep.pit_elev = -std::numeric_limits<double>::infinity();
  }


  //Here we find the pit cells of internally-draining regions. We define these
  //to be cells without any downstream neighbours. Note that this means we will
  //identify all flat cells as being pit cells. For DEMs with many flat cells,
  //this will bloat the priortiy queue slightly. If your DEM includes extensive,
  //predictably located flat regions, you may wish to add these in some special
  //way. Alternatively, you could use Barnes (2014, "An Efficient Assignment of
  //Drainage Direction Over Flat Surfaces") as a way of reducing the number of
  //flat cells. Regardless, the algorithm will deal gracefully with the flats it
  //finds.
  for(int y=1;y<topo.height()-1;y++)
  for(int x=1;x<topo.width() -1;x++){
    const auto my_elev = topo(x,y);
    bool has_lower     = false;
    for(int n=0;n<neighbours;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(topo(nx,ny)<my_elev){
        has_lower = true;
        break;
      }
    }
    //The cell can't drain
    if(!has_lower)
      pq.emplace(x,y,topo(x,y)); 
  }

  //Visit cells in order from lowest to highest
  while(!pq.empty()){
    const auto c = pq.top();                //Copy cell with lowest elevation from priority queue
    pq.pop();                               //Remove the copied cell from the priority queue
    const auto celev = c.z;                 //Elevation of focal cell
    const auto ci    = topo.xyToI(c.x,c.y); //Flat-index of focal cell
    auto clabel      = label(ci);           //Nominal label of cell

    if(clabel==OCEAN){
      //This cell is an ocean cell or a cell that flows into the ocean without
      //encountering any depressions on the way. Upon encountering it we do not
      //need to do anything special.
    } else if(clabel==NO_DEP){
      //Since cells label their neighbours and ocean cells are labeled in the
      //initialization, the only way to get to a cell that is still labeled as
      //not being part of a depression is if that cell were added as a pit cell.
      //For each pit cell we find, we make a new depression and label it
      //accordingly. Not all the pit cells originally added will form new
      //depressions as flat cells will relabel their neighbours.
      clabel          = depressions.size();         //In a 0-based indexing system, size is equal to the id of the next flat
      auto &newdep    = depressions.emplace_back(); //Add the next flat (increases size by 1)
      newdep.pit_cell = topo.xyToI(c.x,c.y);        //Make a note of the pit cell's location
      newdep.pit_elev = celev;                      //Make a note of the pit cell's elevation
      label(ci)       = clabel;                     //Update cell with new label
      // std::cerr<<"\tNew depression from pit cell with label = "<<clabel<<" at "<<c.x<<" "<<c.y<<std::endl;
    } else {
      //Cell has already been assigned to a depression. In this case, it will
      //have no neighbours that can be added to the priority queue so nothing
      //will happen. We use this rather than a boolean array in order to save
      //memory.
    }

    //Consider the cell's neighbours
    for(int n=0;n<neighbours;n++){
      int nx = c.x + mdx[n];              //Get neighbour's coordinates using an offset
      int ny = c.y + mdy[n];              //Get neighbour's coordinates using an offset
      WrapCoordinates(topo,nx,ny);        //Wrap coordinates around torous
      if(!topo.inGrid(nx,ny))             //Is this a valid cell?
        continue;                         //Neighbour cell is out of bounds
      const auto ni     = topo.xyToI(nx,ny);
      const auto nlabel = label(ni);

      if(nlabel==NO_DEP){                 //Neighbour has not been visited yet
        label(ni) = clabel;               //Give the neighbour my label
        pq.emplace(nx,ny,topo(ni));       //Add the neighbour to the priority queue
      } else if (nlabel==clabel) {
        //Skip because we are not interested in ourself. That would be vain.
      } else {
        //We've found a neighbouring depression to which we have not previously
        //found an outlet! Mark the outlet now.

        //Find the elevation and location of the outlet
        auto out_cell = ci;
        auto out_elev = celev;
        if(topo(ni)>out_elev){
          out_cell = ni;
          out_elev = topo(ni);
        }

        outlet_database.emplace(clabel,nlabel,out_cell,out_elev);
      }
    }

  if(label.width()<1000){
    for(int y=0;y<label.height();y++){
      for(int x=0;x<label.width();x++)
        std::cerr<<std::setw(3)<<topo(x,y)<<" ";  
      std::cerr<<"    ";  
      for(int x=0;x<label.width();x++)
        std::cerr<<std::setw(3)<<label(x,y)<<" ";
      std::cerr<<std::endl;
    }
    std::cerr<<std::endl;
  }

  }








  std::vector<Outlet> outlets;

  //Add each outlet. Notice that each outlet is loaded twice. We will fix that
  //below.
  outlets.reserve(outlet_database.size());
  for(const auto &v: outlet_database)
    outlets.push_back(v);

  outlet_database.clear();

  //Sort outlets in order from lowest to highest
  std::sort(outlets.begin(), outlets.end(), [](const Outlet &a, const Outlet &b){
    return a.out_elev<b.out_elev;
  });

  DisjointDenseIntSet djset(depressions.size());
  for(auto &outlet: outlets){
    auto depa_set = djset.findSet(outlet.depa);
    auto depb_set = djset.findSet(outlet.depb);
    // std::cerr<<"Considering "<<outlet.depa<<" "<<outlet.depb<<std::endl;
    // std::cerr<<"\tConsidering "<<depa_set<<" "<<depb_set<<std::endl;
    if(depa_set==depb_set){
      // std::cerr<<"\tSkipping\n";
      continue;
    }
    if(depa_set==OCEAN || depb_set==OCEAN){
      //If we're here then both depressions cannot be the ocean, since we would
      //have used `continue` above. Therefore, one and only one of them is the
      //ocean. We swap them to ensure that the ocean is `depb`
      if(depa_set==OCEAN){
        std::swap(outlet.depa, outlet.depb);
        std::swap(depa_set, depb_set);
      }

      auto &dep = depressions.at(depa_set);

      // std::cerr<<"\tMerging "<<depa_set<<" into the ocean via "<<outlet.depb<<"!"<<std::endl;

      //If this depression has already found the ocean then don't merge it
      //again.
      if(dep.out_cell==OCEAN)
        continue;

      //If this depression already has an outlet, then there's a big problem.
      assert(dep.out_cell==-1);

      //Point this depression to the ocean
      dep.parent   = outlet.depb;
      dep.out_elev = outlet.out_elev;
      dep.out_cell = outlet.out_cell;
      djset.mergeAintoB(depa_set,OCEAN);
    } else {
      //We haven't found the ocean yet, so we merge the two depressions into a
      //new depression.
      auto &depa          = depressions.at(depa_set);
      auto &depb          = depressions.at(depb_set);
      const auto newlabel = depressions.size();
      // std::cerr<<"\tMerging "<<depa_set<<" and "<<depb_set<<" into "<<newlabel<<"!"<<std::endl;
      // std::cerr<<"\tNew parent = "<<newlabel<<std::endl;
      depa.parent         = newlabel;
      depb.parent         = newlabel;
      depa.out_cell       = outlet.out_cell;
      depb.out_cell       = outlet.out_cell;
      depa.out_elev       = outlet.out_elev;
      depb.out_elev       = outlet.out_elev;
      depressions.emplace_back(); //Be sure that this happens AFTER we are done using the `depa` and `depb` references or they may be invalidated!
      djset.mergeAintoB(depa_set, newlabel);
      djset.mergeAintoB(depb_set, newlabel);
    }
  }


  return depressions;
}



void LastLayer(Array2D<label_t> &label, const Array2D<float> &topo, const std::vector<Depression> &depressions){
  #pragma omp parallel for collapse(2)
  for(int y=0;y<label.height();y++)
  for(int x=0;x<label.width();x++){
    auto mylabel = label(x,y);
    while(depressions.at(mylabel).parent!=0 && depressions.at(mylabel).parent!=-1)
      mylabel = depressions.at(mylabel).parent;
    label(x,y) = mylabel;
  }
}



int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  std::string in_name   = argv[1];
  std::string out_name  = argv[2];
  std::string out_graph = argv[3];

  Array2D<float> dem(in_name);   //Recharge (Percipitation minus Evapotranspiration)

  Array2D<label_t> label(dem,NO_DEP);
  const auto deps = GetDepressionHierarchy(dem, label);

  if(label.width()<1000)
  for(int y=0;y<label.height();y++){
    for(int x=0;x<label.width();x++)
      std::cout<<std::setw(3)<<label(x,y)<<" ";
    std::cout<<std::endl;
  }

  std::ofstream fgraph(out_graph);
  fgraph<<"digraph {\n";
  for(unsigned int i=0;i<deps.size();i++)
    fgraph<<i<<" -> "<<deps[i].parent<<";\n";
  fgraph<<"}\n";

  // LastLayer(label, dem, deps);

  label.saveGDAL(out_name);

  return 0;
}