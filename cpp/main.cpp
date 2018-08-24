#include "Array2D.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

//1 2 3
//0   4
//7 6 5
//             0  1  2  3 4 5 6  7              
const int dx8[8] = {-1,-1, 0, 1,1,1,0,-1};
const int dy8[8] = {0, -1,-1,-1,0,1,1, 1};

//  1 
//0   2
//  3
//             0  1 2 3
const int dx4[8] = {-1, 0,1,0};
const int dy4[8] = { 0,-1,0,1};

const int OCEAN_LEVEL = 0;


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
    return z>a.z;
  }
};



double HydroHeight(const int x, const int y, const Array2D<float> &topo, const Array2D<float> &wtd){
  return topo(x,y)+wtd(x,y);
}


void SurfaceWater(const Array2D<float> &topo, Array2D<float> &wtd){
  const int *const dx  = dx8;
  const int *const dy  = dy8;
  const int neighbours = 8;

  Array2D<bool> processed(topo.width(),topo.height(),false);

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
      if(!topo.inGrid(nx,ny)) //Edge cell
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
    std::cout<<"Enqueuing "<<x<<" "<<y<<std::endl;
    q.emplace(x,y,myelev);
  }

  //Move water from high places to low places
  while(!q.empty()){
    //Get the highest cell in the hydroscape
    const auto c = q.top();
    q.pop();

    std::cout<<c.x<<" "<<c.y<<" "<<std::setprecision(40)<<c.z<<std::endl;

    processed(c.x,c.y) = true;

    assert(topo(c.x,c.y)>OCEAN_LEVEL);

    //Figure out where we're sending this cell's water
    double max_slope = 0;   //Maximum slope we've seen so far. Setting to 0 ensures that we only consider downhill neighbours.
    int max_n        = -1;  //Lowest downhill neighbour
    for(int n=0;n<neighbours;n++){
      const int nx       = c.x + dx[n];
      const int ny       = c.y + dy[n];
      const auto nheight = HydroHeight(nx,ny,topo,wtd);
      if(!topo.inGrid(nx,ny)) //Edge cell
        continue;

      const double slope = c.z-nheight; //Slope to the neighbour
      if(slope>max_slope){
        max_slope = slope;
        max_n     = n;
      }
    }

    //max_n is the lowest neighbour
    if(max_n==-1) //There was no lowest neighbour
      continue;

    const int  nx      = c.x+dx[max_n];
    const int  ny      = c.y+dy[max_n];
    const auto nheight = HydroHeight(nx,ny,topo,wtd);

    std::cout<<"\tLowest neighbour: "<<nx<<" "<<ny<<std::endl;

    //We've found a downhill neighbour, now we need to move water that direction

    //Is my neighbour an ocean?
    if(topo(nx,ny)<=OCEAN_LEVEL){ //TODO: This might be a special value
      wtd(c.x,c.y) = 0;
      wtd(nx,ny)   = 0; //TODO: Probably unnecessary
    } else {
      const double water_to_move = std::min(wtd(c.x,c.y),(float)((c.z-nheight)/2.0));
      if(water_to_move>1e-4){
        wtd(c.x,c.y) -= water_to_move;
        wtd(nx,ny)   += water_to_move;   
        std::cout<<"\tWater to move: "<<water_to_move<<std::endl;

        //Enqueue the neighbour we've just passed water to
        std::cout<<"Enqueuing "<<nx<<" "<<ny<<std::endl;
        q.emplace(nx,ny,HydroHeight(nx,ny,topo,wtd));
      }
    }

    const auto cheight_new = HydroHeight(c.x,c.y,topo,wtd);

    //Add appropriate neighbours to priority queue
    for(int n=0;n<neighbours;n++){
      const int  nx      = c.x + dx[n];
      const int  ny      = c.y + dy[n];
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
      if(processed(nx,ny)) // && std::abs(nheight-cheight_new)<1e-6)
        continue;

      std::cout<<"Enqueuing "<<nx<<" "<<ny<<std::endl;
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
}