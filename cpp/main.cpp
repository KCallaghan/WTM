#include "Array2D.hpp"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

//1 2 3
//0   4
//7 6 5
//                      0  1  2  3 4 5 6  7              
const int dx8[8]    = {-1,-1, 0, 1,1,1,0,-1};
const int dy8[8]    = {0, -1,-1,-1,0,1,1, 1};
const int d8inverse = {4,  5, 6, 7,0,1,2, 3};

//  1 
//0   2
//  3
//                      0  1 2 3
const int dx4[8]    = {-1, 0,1,0};
const int dy4[8]    = { 0,-1,0,1};
const int d4inverse = { 2, 3,0,1};

const int *const dx       = dx8;
const int *const dy       = dy8;
const int *const dinverse = d8inverse;
const int neighbours      = 8;


const int OCEAN_LEVEL = 0;

class GridCell {
 public:
  int x, y;
  GridCell(const int x0, const int y0){
    x = x0;
    y = y0;
  }
};


class GridCellZ {
 public:
  int x, y;
  float z;
  GridCellZ(const int x0, const int y0, const float z0){
    x = x0;
    y = y0;
    z = z0;
  }
  bool operator>(const GridCellZ& a) const { 
    return z>a.z; //Less than sorts the queue in reverse
  }
};



double HydroHeight(const int x, const int y, const Array2D<float> &topo, const Array2D<float> &wtd){
  return topo(x,y)+wtd(x,y);
}



int GetLowestNeighbour(const int x, const int y, const Array2D<float> &topo, const Array2D<float> &wtd){
  //Figure out where we're sending this cell's water
  double max_slope = 0;   //Maximum slope we've seen so far. Setting to 0 ensures that we only consider downhill neighbours.
  int    max_n     = -1;  //Lowest downhill neighbour
  const auto myheight = HydroHeight(x,y,topo,wtd);
  for(int n=0;n<neighbours;n++){
    const int nx       = x + dx[n];
    const int ny       = y + dy[n];
    const auto nheight = HydroHeight(nx,ny,topo,wtd);
    if(!topo.inGrid(nx,ny)) //Edge cell
      continue;

    const double slope = myheight-nheight; //Slope to the neighbour
    if(slope>max_slope){
      max_slope = slope;
      max_n     = n;
    }
  }

  return max_n;
}



void ComplicatedDepressionDistribute(){

}


//Returns the amount of excess water after filling the depression. May be 0.
double DepressionDistribute(
  const int deplabel,
  const Array2D<int> &label,
  const double spill_elevation,
  Array2D<float> &wtd,
  double water_volume
){
  //Find volume of depression
  double dep_volume = 0;
  #pragma omp parallel for collapse(2) reduction(+:dep_volume)
  for(int y=0;y<wtd.height();y++)
  for(int x=0;x<wtd.width();x++){
    if(label(x,y)!=deplabel)
      continue;
    dep_volume += spill_elevation - topo(x,y);
    if(wtd(x,y)<0)
      dep_volume += -wtd(x,y);
  }

  if(dep_volume<water_volume){
    //Depression couldn't hold all the incoming water, so we fill it entirely.
    #pragma omp parallel for collapse(2)
    for(int y=0;y<wtd.height();y++)
    for(int x=0;x<wtd.width();x++){   
      if(label(x,y)!=deplabel)
        continue;
      wtd(x,y) = spill_elevation - topo(x,y);
    } 
    return water_volume-dep_volume;
  } else {
    ComplicatedDepressionDistribute()
    return 0; //TODO
  }
}


void SurfaceWater2(const Array2D<float> &topo, Array2D<float> &wtd){
  const int NOT_DEP = -1;
  Array2D<int>  label    (topo.width(),topo.height(),NOT_DEP); //All depressions are unlabeled
  Array2D<bool> processed(topo.width(),topo.height(),false);
  Array2D<int>  flowdirs (topo.width(),topo.height(),-1   );

  int next_dep_id = 0;                //First value that is not unlabaled

  std::vector<float> spill_elevation;
  std::map<GridCell, double> depression_volume;

  //Sort cells so that the lowest cell comes off the queue first
  std::priority_queue<GridCellZ, std::vector<GridCellZ>, std::greater<GridCellZ> > pq;

  //Add all edge cells to the queue TODO: Terribly inefficient
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++)
    if(topo.isEdgeCell(x,y))
      pq.emplace(x,y,topo(x,y));


  while(!pq.empty()){
    const auto c = pq.top();
    pq.pop();

    for(int n=0;n<neighbours;n++){
      const int nx = x + dx[n];
      const int ny = y + dy[n];
      if(!topo.inGrid(nx,ny)) //Cell is out of bounds
        continue;
      if(processed(nx,ny))    //Cell has already been in the queue
        continue;
      processed(nx,ny) = true;

      //Is the neighbouring cell lower than I am?
      if(topo(nx,ny)<=topo(x,y)){ //Cell is part of a depression
        int deplabel = label(x,y);
        if(deplabel==0){ //I don't have a label: this is a new depression
          deplabel = next_dep_id++;
          spill_elevation.push_back(topo(x,y));
          depression_volume.push_back(0);
        }
        label(nx,ny) = deplabel;
        pq.emplace(nx,ny,topo(x,y));
      } else {
        pq.emplace(nx,ny,topo(nx,ny));
      }
      flowdirs(nx,ny) = dinverse[n];
    }
  }

  //At this point every depression is labeled with a unique id [0,*), each id
  //has a corresponding spill elevation, and each cell has a pointer to its
  //downstream "parent" cell.

  //Calculate how many upstream cells flow into each cell
  Array2D<char>  dependencies(topo.width(),topo.height(),0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++)
  for(int n=0;n<neighbours;n++)
    if(flowdirs[n]==dinverse[n])
      dependencies(x,y)++;

  //Find the peaks
  std::queue<GridCell> q;
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++)
    if(dependencies(x,y)==0)
      q.emplace(x,y);

  while(!q.empty()){
    const auto c = q.front();
    q.pop();

    const auto n = flowdirs(c.x,c.y);
    const int  nx = c.x+dx[n];
    const int  ny = c.y+dy[n];

    if(label(c.x,c.y)==NOT_DEP && label(nx,ny)==NOT_DEP){ //I'm not a depression and neither is downstream neighbour
      if(wtd(c.x,c.y)>=0){
        wtd(nx,ny)  += wtd(c.x,c.y);
        wtd(c.x,c.y) = 0;
      }
    } else if(label(c.x,c.y)==NOT_DEP && label(nx,ny)!=NOT_DEP){ //Passing flow into a depression
      if(wtd(c.x,c.y)>=0){
        depression_volume[label(nx,ny)] += wtd(c.x,c.y);
        wtd(c.x,c.y)                     = 0;
      }
    } else if(label(c.x,c.y)!=NOT_DEP && label(nx,ny)==NOT_DEP){ //Passing flow out of a depression
      const auto this_dep        = label(c.x,c.y);
      const auto leftover_volume = DepressionDistribute(
        this_dep,
        label,
        spill_elevation[this_dep],
        wtd,
        depression_volume[this_dep]
      );
      wtd(nx,ny) += leftover_volume;
    } else { //Passing flow around a depression
      if(topo(nx,ny)>topo(c.x,c.y)){ //This is the pit cell of a depression
        depression_volume[c] = wtd(c.x,c.y);
        wtd(c.x,c.y) = 0;
      } else {                       //Moving towards a pit cell
        if(wtd(c.x,c.y)>=0){
          wtd(nx,ny)  += wtd(c.x,c.y);
          wtd(c.x,c.y) = 0;
        }
      }
    }

    if(--dependencies(nx,ny)==0)
      q.emplace(nx,ny);
  }
}




void SurfaceWater(const Array2D<float> &topo, Array2D<float> &wtd){
  Array2D<int> processed(topo.width(),topo.height(),0);

  std::priority_queue<GridCellZ, std::vector<GridCellZ>, std::greater<GridCellZ> > q;

  //Find the hydrologic peaks: those places into which no 
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    bool has_higher = false;
    const auto myelev = HydroHeight(x,y,topo,wtd); //Elevation of focal cell

    //The current cell has no water, so we ignore it.
    if(wtd(x,y)<=0)
      continue;
    //The current cell is an ocean, so ignore it.
    if(topo(x,y)<=OCEAN_LEVEL)
      continue;

    //Is any neighbour higher than me?
    for(int n=0;n<neighbours;n++){
      const int nx = x + dx[n];
      const int ny = y + dy[n];
      if(!topo.inGrid(nx,ny)) //Cell is out of bounds
        continue;
      //Neighbour cell has no water, so ignore it. It may be higher than the
      //focal cell, but it can't pass flow into the focal cell.
      if(wtd(nx,ny)<=0)             
        continue;
      const auto nelev = HydroHeight(nx,ny,topo,wtd);
      if(nelev>myelev){
        has_higher = true;
        break;
      }
    }

    //I had a higher neighbour, so I am not a peak
    if(has_higher)
      continue;

    //I am a peak
    std::cerr<<"Enqueuing "<<x<<" "<<y<<std::endl;
    q.emplace(x,y,myelev);
    processed(x,y)++;
  }

  std::cout<<"Peaks found = "<<q.size()<<std::endl;

  //Move water from high places to low places
  while(!q.empty()){
    //Get the highest cell in the hydroscape
    const auto c = q.top();
    q.pop();

    //Cell has been modified since it was emplaced, so we skip it
    if(HydroHeight(c.x,c.y,topo,wtd)!=c.z)
      continue;

    std::cerr<<"Popping "<<c.x<<" "<<c.y<<" "<<std::setprecision(40)<<c.z<<std::endl;
    // std::cout<<"\tNew top: "<<q.top().x<<" "<<q.top().y<<std::endl;

    processed(c.x,c.y)++;

    // if(processed(c.x,c.y)>100)
      // continue;

    assert(topo(c.x,c.y)>OCEAN_LEVEL);

    const auto max_n = GetLowestNeighbour(c.x,c.y,topo,wtd);

    //max_n is the lowest neighbour
    if(max_n==-1) //There was no lowest neighbour
      continue;

    const int  nx      = c.x+dx[max_n];
    const int  ny      = c.y+dy[max_n];
    const auto nheight = HydroHeight(nx,ny,topo,wtd);

    // std::cerr<<"\tLowest neighbour: "<<nx<<" "<<ny<<std::endl;

    //We've found a downhill neighbour, now we need to move water that direction

    //Is my neighbour an ocean?
    if(topo(nx,ny)<=OCEAN_LEVEL){ //TODO: This might be a special value
      wtd(c.x,c.y) = 0;
      wtd(nx,ny)   = 0; //TODO: Probably unnecessary
    } else {
      const double water_to_move = std::min(wtd(c.x,c.y),(float)((c.z-nheight)/2.0));
      if(water_to_move>1e-3){
        wtd(c.x,c.y) -= water_to_move;
        wtd(nx,ny)   += water_to_move;
        std::cerr<<"\tWater to move: "<<std::setprecision(40)<<water_to_move<<std::endl;

        //Enqueue the neighbour we've just passed water to
        std::cerr<<"\tEnqueuing (move) "<<nx<<" "<<ny<<std::endl;
        q.emplace(nx,ny,HydroHeight(nx,ny,topo,wtd));
      }
    }

    const auto cheight_new = HydroHeight(c.x,c.y,topo,wtd);

    //Add appropriate neighbours to priority queue
    for(int n=0;n<neighbours;n++){
      const int  nx      = c.x + dx[n];
      const int  ny      = c.y + dy[n];
      // std::cout<<"\tConsidering "<<nx<<" "<<ny<<std::endl;
      if(!topo.inGrid(nx,ny))
        continue;

      const auto nheight = HydroHeight(nx,ny,topo,wtd);

      //Don't add the neighbour we passed water to, since we added it already
      //above
      if(n==max_n)
        continue;
      //If there's no water in the neighbour, we can skip it.
      if(wtd(nx,ny)<=0)
        continue;
      //Don't process the oceans
      if(topo(nx,ny)<=OCEAN_LEVEL)
        continue;
      //Don't add neighbours whose hydrologic height is too similar to the focal
      //cell
      if(processed(nx,ny)>0) // && std::abs(nheight-cheight_new)<1e-6)
        continue;

      std::cerr<<"Enqueuing (neigh) "<<nx<<" "<<ny<<std::endl;
      processed(nx,ny)++;
      q.emplace(nx,ny,nheight);
    }
  }
}





int main(int argc, char **argv){
  if(argc!=3){
    std::cout<<"Syntax: "<<argv[0]<<" <Input file directory> <Run Type>"<<std::endl;
    return -1;
  }

  std::string dir      = argv[1];
  std::string run_type = argv[2];

  Array2D<float> rech  (dir+"/Mad_020500_rech_rotated.nc",   "value");   //Recharge (Percipitation minus Evapotranspiration)
  Array2D<float> temp  (dir+"/Mad_020500_temp_rotated.nc",   "value");   //Air temperature   - Used with fslope to generate efolding depth
  Array2D<float> fslope(dir+"/Mad_020500_fslope_rotated.nc", "value");   //100/(1+150*slope) - Used with fslope to generate efolding depth
  Array2D<float> topo  (dir+"/Mad_020500_topo_rotated.nc",   "value");   //Terrain height
  Array2D<float> ksat  (dir+"/Mad_ksat_rotated.nc",          "value");   //Hydrologic conductivity
  Array2D<float> wtd   (dir+"/Mad_021000_wtd_rotated.nc",    "value"); 

  if(run_type=="equilibrium")
    throw std::runtime_error("equilibrium not implemented!");
  else if (run_type=="transient"){
    //Pass
  } else 
    throw std::runtime_error("Expected 'equilibrium' or 'transient'!");

  //TODO: Make sure all files have same dimensions

  for(int i=0;i<wtd.size();i++){
    wtd(i) = 0;
    if(topo(i)>0)
      wtd(i) = 1;
  }

  SurfaceWater(topo, wtd);

  std::ofstream fout("/z/out.dat");
  fout.write(reinterpret_cast<const char*>(wtd.data), wtd.size()*sizeof(float));

  return 0;
}