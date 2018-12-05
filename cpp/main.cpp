//I can't figure out how to use the test cases (throws errors when trying to open the file)

#include "dephier_b.hpp"
#include "DisjointDenseIntSet.hpp"
#include "netcdf.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <richdem/common/Array2D.hpp>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace rd = richdem;

const int    *const dx       = dx8;
const int    *const dy       = dy8;
const int    *const dinverse = d8inverse;
const double *const dr       = dr8;
const int neighbours         = 8;


const float  OCEAN_LEVEL = -9999;  //ocean_level in the topo file must be lower than any non-ocean cell. 

rd::Array2D<flowdir_t> flowdirs; //TODO: Make non-global
template<class T>
void PrintStamp(const std::string title, const rd::Array2D<T> &arr, const int x0, const int y0, const int r, const int width=2){
  std::cout<<"\n"<<title<<std::endl;
  for(int y=std::max(y0-r,0);y<std::min(y0+r+1,arr.height());  y++){
    for(int x=std::max(x0-r,0);x<std::min(x0+r+1,arr.width()); x++){
      if(x==x0 && y==y0) std::cout<<"\033[91m";
      std::cout<<std::setw(width)<<arr(x,y)<<" ";
      if(x==x0 && y==y0) std::cout<<"\033[39m";
    }
    std::cout<<std::endl;
  }
}

template<class T>
void PrintStampWithSpecialCells(const std::string title, const rd::Array2D<T> &arr, const std::vector<int> &sc, const int x0, const int y0, const int r, const int width=2){
  std::cout<<"\n"<<title<<std::endl;
  for(int y=std::max(y0-r,0);y<std::min(y0+r+1,arr.height());  y++){
    for(int x=std::max(x0-r,0);x<std::min(x0+r+1,arr.width()); x++){
      bool colour = false;
      if(x==x0 && y==y0)                   {std::cout<<"\033[91m"; colour=true;}
      else if(std::find(sc.begin(),sc.end(),arr.xyToI(x,y))!=sc.end()) {std::cout<<"\033[92m"; colour=true;}
      std::cout<<std::setw(width)<<arr(x,y)<<" ";
      if(colour) std::cout<<"\033[39m";
    }
    std::cout<<std::endl;
  }
}


template<class T>
void PrintDEM(const std::string title, const rd::Array2D<T> &arr, const int width=2){
  std::cout<<"\n"<<title<<std::endl;
  for(int y=0;y<arr.height();y++){
    for(int x=0;x<arr.width(); x++)
      std::cout<<std::setw(width)<<arr(x,y)<<" ";
    std::cout<<std::endl;
  }
}

template<>
void PrintDEM(const std::string title, const rd::Array2D<flowdir_t> &arr, const int width){
  std::cout<<"\n"<<title<<std::endl;
  for(int y=0;y<arr.height();y++){
    for(int x=0;x<arr.width(); x++)
      std::cout<<std::setw(width)<<(int)arr(x,y)<<" ";
    std::cout<<std::endl;
  }
}


//TODO
// template<class T>
// void mWvol(const int dep, std::vector<Depression<elev_t> > &deps){
// 
// }


//Richard: Checked this
template<class elev_t>
void SurfaceWater(
  const rd::Array2D<elev_t>        &topo,
  rd::Array2D<float>               &wtd,
  const rd::Array2D<int>           &label,
  std::vector<Depression<elev_t> > &deps,
  const rd::Array2D<flowdir_t>     &flowdirs
){
  //Our first step is to move all of the water downstream into pit cells. To do
  //so, we use the steepest-descent flow directions provided by the depression
  //hierarchy code

  //Calculate how many upstream cells flow into each cell
  rd::Array2D<char>  dependencies(topo.width(),topo.height(),0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width(); x++){
  for(int n=0;n<neighbours;n++){                           //Loop through neighbours
    const int nx = x+dx[n];                                //Identify coordinates of neighbour
    const int ny = y+dy[n];
    if(!topo.inGrid(nx,ny))
      continue;    
    if(flowdirs(nx,ny)==dinverse[n])  //CHECK WITH RICHARD           //(int)topo.xyToI(x,y))              //Does my neighbour flow into me?
      dependencies(x,y)++;                                 //Increment my dependencies
  }
  //std::cout<<"cell "<<x<<" "<<y<<" dependencies "<<dependencies(x,y)<<std::endl;
}

  int pit_cell_count = 0;
  int peak_count     = 0;
  int flat_count     = 0;
  for(unsigned int i=0;i<flowdirs.size();i++){
    if(dependencies(i)!=0 && flowdirs(i)==NO_FLOW)
      pit_cell_count++;
    if(dependencies(i)==0 && flowdirs(i)!=NO_FLOW)
      peak_count++;
    if(dependencies(i)==0 && flowdirs(i)==NO_FLOW)
      flat_count++;
  }
  std::cout<<"Found "<<pit_cell_count<<" pit cells."<<std::endl;
  std::cout<<"Found "<<peak_count    <<" peak cells."<<std::endl;
  std::cout<<"Found "<<flat_count    <<" flat cells."<<std::endl;


  //Find the peaks. These are the cells into which no other cells pass flow (i.e. 0 dependencies). We
  //know the flow accumulation of the peaks without having to perform any
  //recursive calculations; they just pass flow downstream. From the peaks, we
  //can begin a breadth-first traversal in the downstream direction by adding
  //each cell to the frontier/queue as its dependency count drops to zero.
  std::queue<int> q;
  for(unsigned int i=0;i<topo.size();i++){
    if(dependencies(i)==0)// && flowdirs(i)!=NO_FLOW)  //Is it a peak?
      q.emplace(i);       
  }  //Yes.

  int pit_cells_found = 0; //TODO: For debugging
  int cells_traversed = 0; //TODO: For debugging

  //Starting with the peaks, pass flow downstream
  while(!q.empty()){
    const auto c = q.front();          //Copy focal cell from queue
    q.pop();                           //Clear focal cell from queue

    cells_traversed++; //TODO: For debugging

    //Coordinates of downstream neighbour, if any
    const auto ndir = flowdirs(c); 

    int n = NO_FLOW;
    if(ndir==NO_FLOW){ //TODO: For debugging
      pit_cells_found++;
    } else { //TODO: Fix this monkey patching
      int x,y;
      topo.iToxy(c,x,y);
      const int nx = x+dx[ndir];
      const int ny = y+dy[ndir];
      n            = topo.xyToI(nx,ny);
      assert(n>=0);
    }

    //TODO: Might need this - could also check label
    //If downstream neighbour is the ocean, we drop our water into it and the
    //ocean is unaffected.
    // if(wtd(c)<=OCEAN_LEVEL){                                                    //I am confused how this is actually checking if downstream neighbour is the ocean. OCEAN_LEVEL = 0, so it looks more like here is a place to reset any accidentally negative wtd to 0? When coupled, negative wtd would actually be allowed, so we shouldn't do this
    //   wtd(c) = 0;
    //   continue;
    // }
  
    if (n == NO_FLOW){    //if this is a pit cell, move the water to the appropriate depression's water_vol.    
      if(wtd(c)>0){
        deps[label(c)].water_vol += wtd(c);
        // std::cout<<"Giving water to "<<label(c)<<" "<<c<<" "<<deps[label(c)].water_vol<<std::endl;
   //   std::cout<<"the total water in this depression is "<<deps[label(c)].water_vol<<" and it is depression "<<label(c)<<std::endl;
        wtd(c) = 0; //Clean up as we go
      }
    } else {                               //not a pit cell
      //If we have water, pass it downstream.
      if(wtd(c)>0){ //Groundwater can go negative, so it's important to make sure that we are only passing positive water around
    //    std::cout<<"we have water coming from "<<c<<" and going to "<<n<<std::endl;
        wtd(n) += wtd(c);  //Add water to downstream neighbour. This might result in filling up the groundwater table, which could have been negative
        wtd(c)  = 0;       //Clean up as we go
      }
  
      //Decrement the neighbour's dependencies. If there are no more dependencies,
      //we can process the neighbour.
      if(--dependencies(n)==0){                            //CHECK this is a hack! Something is wrong with the dependencies matrix, should never go below 0 but it sometimes does. 
        assert(dependencies(n)>=0);
        q.emplace(n);                   //Add neighbour to the queue
      }
    }
  }

  std::cout<<"Found pit cells = "<<pit_cells_found<<std::endl;
  std::cout<<"Cells traversed = "<<cells_traversed<<std::endl;
}


 //TODO: 0. Calculate the number of cells within and volume of each depression.   DONE 
  //This will take place inside of GetDepressionHierarchy.
 

 //TODO 1. DONE. Get flow directions for all cells. 

  //TODO 2. Perform a flow accumulation moving water downhill and filling             NEARLY DONE - think about oceans
  //groundwater as you go. Look at `misc/main.cpp`. Recall that we did this by
  //counting depressions, finding peaks, and using a queue to control a breadth-
  //first traversal from the peaks downwards.



  //TODO 3. As part of the above, when flow can't go downhill any farther, add        DONE
  //it to the `water_vol` for the appropriate depression. Use the labels array
  //to determine the appropriate depression.




//Richard: Checked this
template<class elev_t>
void Overflow(
  int current_depression,
  std::vector<Depression<elev_t> > &deps
){
  if(current_depression==NO_VALUE)
    return;

  auto &this_dep = deps.at(current_depression);

  //Visit child depressions. When these both overflow, then we spread water
  //across them by spreading water across their common metadepression
  Overflow(this_dep.lchild, deps);
  Overflow(this_dep.rchild, deps);

  //Catch depression that link to the ocean through this one. These are special
  //cases because we will never spread water across the union of these
  //depressions and the current depressions: they only flow into the current
  //depression
  for(const auto c: this_dep.ocean_linked)
    Overflow(c, deps);

  //If the current depression is the ocean then at this point we've visited all
  //of its ocean-linked depressions (the ocean has no children). Since we do not
  //otherwise want to modify the ocean we quit here.
  if(current_depression==OCEAN)
    return;

  // std::cout<<"depression number "<<this_dep.dep_label<<" volume "<<this_dep.dep_vol<<" water "<<this_dep.water_vol<<std::endl;

  //Each depression has an associated dep_vol. This is the TOTAL volume of the
  //meta-depression including all of its children. This property answers the
  //question, "How much water can this meta-depression hold?"

  //Each depression also has an asoociated water volume called `water_vol`. This
  //will stored the MARGINAL volume of the water in the depression. That is, for
  //each depression this indicates how much water is in this depression that
  //would not fit into its child depressions.

  //We are overflowing the depression
  if(this_dep.odep == OCEAN){
    //If a depression overflows directly into an ocean then its odep is the
    //ocean and so is its parent.

    //The current depression's outlet is into the ocean. Since the ocean can
    //absorb an infinite amount of water without changing its water volume, we
    //simply set the amount of water contained in the current depression to be
    //either its water_vol (the depression doesn't overflow) or equal to its
    //depression volume (the excess water is thrown into the ocean, which is
    //unaffected).
    this_dep.water_vol = std::min(this_dep.water_vol,this_dep.dep_vol);
  } else if(this_dep.water_vol>this_dep.dep_vol) {
    //The neighbouring depression is not the ocean and this depression is
    //overflowing (therefore, all of its children are full)
    assert(this_dep.lchild==NO_VALUE || deps.at(this_dep.lchild).water_vol==deps.at(this_dep.lchild).dep_vol);
    assert(this_dep.rchild==NO_VALUE || deps.at(this_dep.rchild).water_vol==deps.at(this_dep.rchild).dep_vol);

    //The marginal volume of this depression is larger than what it can hold, so
    //we determine the amount that overflows, the "extra water".
    double extra_water  = this_dep.water_vol - this_dep.dep_vol;
  
    //Now that we've figured out how much excess water there is, we fill this
    //depression to its brim. Note that we don't use addition or subtraction
    //here to ensure floating-point equality.
    this_dep.water_vol = this_dep.dep_vol;

    auto &outlet_dep = deps.at(this_dep.odep);                             //Depression this one overflows into

    outlet_dep.water_vol += extra_water;                                   //Add excess water to the overflow depression
    if(outlet_dep.water_vol>outlet_dep.dep_vol){                           //Water in overflow depression exceeds its volume
      extra_water          = outlet_dep.water_vol - outlet_dep.dep_vol;    //Find out how much overflow depression overflowed by
      outlet_dep.water_vol = outlet_dep.dep_vol;                           //Mark overflow depression as being entirely filled
    }

    //Extra water should always be >=0. If, for floating-point reasons it is a
    //small negative number, we will cross that bridge when we come to it.
    assert(extra_water>=0);

    //Since if a depression links directly to an ocean its odep and parent
    //values are both OCEAN, we have overflowed into the ocean in the upper part
    //of this if-clause.
    assert(this_dep.parent!=OCEAN);

    //At this point extra_water is the amount of water that is left when both
    //this depression and its overflow depression are completely filled. We add
    //this extra water to this depression's (and therefore it's overflow
    //depression's) parent.
    deps.at(this_dep.parent).water_vol += extra_water;
    // std::cout<<"parent number "<<deps.at(this_dep.parent).dep_label<<" now has water volume "<<deps.at(this_dep.parent).water_vol<<std::endl;
  }

  // std::cout<<"and after: depression number "<<this_dep.dep_label<<" volume "<<this_dep.dep_vol<<" water "<<this_dep.water_vol<<std::endl;

  //All overflowing depressions should by now have overflowed all the way down
  //to the ocean. We must now spread the water in the depressions out by setting
  //appropriate values for wtd.
}

 





//Simple data structure to hold information needed to spread water in a filled
//metadepression.
class SubtreeDepressionInfo {
 public:
  //One of the depressions at the bottom of the meta-depression. We use this to
  //identify a pit cell from which to start flooding.
  int   leaf_label = -1;          
  //The metadepression containing all of the children. This metadepression is
  //guaranteed to be large enough to hold all of the water of its children plus
  //whatever exists only in the metadepression itself. We use this to determine
  //the water and depression volumes.
  int   top_label = -1;
  //Here we keep track of which depressions are contained within the
  //metadepression. This allows us to limit the spreading function to cells
  //within the metadepression.
  std::unordered_set<int> my_labels;
  //Total water of all the depressions up to and including the one identified by
  //top_label
  double water_vol = 0;
};




template<class elev_t>
void Fill_Water(
  //Identifies a meta-depression through which water should be spread, leaf node
  //from which the water should be spread, valid depressions across which water
  //can spread, and the amount of water to spread
  SubtreeDepressionInfo                  &stdi,  
  const std::vector<Depression<elev_t> > &deps,  //Depression hierarchy
  const rd::Array2D<float>               &topo,  //Topographic data for calculating marginal volumes as we attempt to spread water
  const rd::Array2D<label_t>             &label, //2D array in which each cell is labeled with the leaf depression it belongs to
  rd::Array2D<float>                     &wtd    //Water table depth: we transfer water into this
){
  //Nothing to do if we have no water
  if(stdi.water_vol==0)
    return;

  //changing tactics to start always from the leaves, then work your way up until you find something that isn't completely full. 
  std::cerr<<"\n\n\033[35m####################### Fill Water\033[39m"<<std::endl;

  //TODO: Use hashset to avoid allocating massive chunks of memory.
  rd::Array2D<bool> visited(topo.width(),topo.height(),false);
 
  //Priority queue that sorts cells by lowest elevation first. If two cells are
  //of equal elevation the one added most recently is popped first. The ordering
  //of the cells processed by this priority queue need not match the ordering of
  //the cells processed by the depression hierarchy.
  GridCellZk_high_pq<elev_t> flood_q;                          

  std::cerr<<"Bottom label      = "<<stdi.leaf_label<<std::endl;
  std::cerr<<"Depression volume = "<<deps.at(stdi.top_label).dep_vol<<"\n";
  std::cerr<<"Water volume      = "<<stdi.water_vol<<std::endl;
  std::cerr<<"Allowed labels    = ";
  for(auto x:stdi.my_labels)
    std::cerr<<x<<" ";
  std::cerr<<std::endl;

 
  { //Scope to limit pit_cell
    //Cell from which we can begin flooding the meta-depression. Which one we
    //choose is arbitrary, since we will fill all of the leaf depressions and
    //intermediate metadepressions until and including when we reach the
    //depression identified by stdi.top_label
    const auto pit_cell    = deps.at(stdi.leaf_label).pit_cell;
    std::cout<<"pit cell "<<pit_cell<<" "<<stdi.leaf_label<<std::endl;
    assert(pit_cell>=0);

    flood_q.emplace(
      pit_cell % topo.width(),
      pit_cell / topo.width(),
      topo(pit_cell)
    );                    //create a new priority queue starting with the pit cell of the depression

    visited(pit_cell) = true;//label(pit_cell);         // show that we have already added this cell to those that have water. We need a better way to do this. 
  }

  double current_volume;              //TODO: This is out here for debugging purposes. SHould be moved inside pq
  GridCellZk_high<elev_t> c(0,0,0,0); //TODO: Out for debugging should be in pq

  //Cells whose wtd will be affected as we spread water around
  std::vector<int> cells_affected;

  //Stores the sum of the elevations of all of the cells in cells_affected. Used
  //for calculating the volume we've seen so far. (See explanation above or in
  //dephier.hpp TODO)
  double total_elevation = 0;


  //TODO: It's possible to accelerate this by greedily eating cells which belong
  //"only" to the meta-depression.

  std::cout<<"about to start the flood_q while loop "<<label(1,1)<<std::endl;
  while(!flood_q.empty()){
    c = flood_q.top(); //TODO local var
    flood_q.pop();

    std::cout<<"width "  <<topo.width()<<std::endl;
    std::cout<<"x and y "<<c.x<< " "<<c.y<<std::endl;
    std::cout<<"current cell elevation "<<topo(c.x,c.y)<<std::endl;

    //We keep track of the current volume of the depression by noting the total
    //elevation of the cells we've seen as well as the number of cells we've
    //seen.

    //TODO: Note that the current cell's above ground volume and wtd do not
    //contribute at all. This a choice that Kerry and Richard discussed. It is
    //as though there is a virtual water line coincident with the edge of the
    //current cell. No water infiltrates into this cell or is stored above it -
    //only cells previously visited are considered when doing volume
    //calculations.

    //Current volume of this subset of the metadepression. Since we might climb
    //over a saddle point, this value can occasionally be negative. It will be
    //positive by the time we need to spread the water around.
    current_volume = cells_affected.size()*topo(c.x,c.y) - total_elevation; //TODO: Local var

    std::cout<<"volume "<<current_volume<<" cells "<<cells_affected.size()<<" topo "<<topo(c.x,c.y)<<" total "<<total_elevation<<std::endl;

    //TODO: If this is false by a small margin, then it's a floating point issue
    //and this should be adjusted to be >=-1e-6 and water_vol should be made 0
    assert(stdi.water_vol>=0); 

    //All the cells within this depression should have water table depths less
    //than or equal to zero because we have moved all of their water down slope
    //into the pit cell. Since we may already have filled other depressions
    //their cells are allowed to have wtd>0. Thus, we raise a warning if we are
    //looking at a cell in this unfilled depression with wtd>0.
    if(stdi.my_labels.count(label(c.x,c.y))==1 && wtd(c.x,c.y)>0){
      PrintDEM("Flowdirs", flowdirs, 9);
      PrintDEM("wtd", wtd, 9);
      PrintDEM("Labels", label, 9);
      throw std::runtime_error("A cell was discovered in an unfilled depression with wtd>0!");
    }

    //There are two possibilities:
    //1. The virtual water level exceeds slightly the height of the cell. The cell's water table then fills up as much as it can.
    //   The water surface is then level with the height of the cell.
    //
    //2. The water surface is below the height of the cell because there is sufficient topographic volume to hold all the water.
    //   In this case, the cell's water table is left unaffected.

    if(stdi.water_vol<=current_volume-wtd(c.x,c.y)){
      //The current scope of the depression plus the water storage capacity of
      //this cell is sufficient to store all of the water. We'll stop adding
      //cells and fill things now.
      const auto my_elev = topo(c.x,c.y);

      std::cerr<<"Attempting to fill depression..."<<std::endl;
      std::cerr<<"\tLabel of last cell       = "<<label(c.x,c.y)       <<std::endl;
      std::cerr<<"\tWater volume             = "<<stdi.water_vol       <<std::endl;
      std::cerr<<"\tDepression volume        = "<<deps.at(stdi.top_label).dep_vol         <<std::endl;
      std::cout<<"\tDepression number        = "<<stdi.leaf_label       <<std::endl;
      std::cerr<<"\tCurrent volume           = "<<current_volume       <<std::endl;
      std::cerr<<"\tTotal elevation          = "<<total_elevation      <<std::endl;
      std::cerr<<"\tCurrent elevation        = "<<my_elev              <<std::endl;
      std::cerr<<"\tNumber of cells affected = "<<cells_affected.size()<<std::endl;  

      //We will fill the depression so that the surface of the water is at this
      //elevation.
      double water_level;

      if(current_volume<stdi.water_vol){ //TODO: Check stdi.my_labels.count(label(c.x,c.y))==0 ?
        //The volume of water exceeds what we can hold above ground, so we will
        //stash as much as we can in this cell's water table. This is okay
        //because the above ground volume plus this cell's water table IS enough
        //volume (per the if-clause above).
        std::cout<<"we are in the if"<<std::endl;
        //Fill in as much of this cell's water table as we can
        const double fill_amount = stdi.water_vol - current_volume;
        assert(fill_amount>=0);
        wtd(c.x,c.y)   += fill_amount;
        stdi.water_vol -= fill_amount;   //Doesn't matter because we don't use water_vol anymore
        water_level     = topo(c.x,c.y);
      } else if (current_volume==stdi.water_vol){
        //The volume of water is exactly equal to the above ground volume so we
        //set the water level equal to this cell's elevation
        water_level = topo(c.x,c.y);
      } else {
        //The water volume is less than this cell's elevation, so we calculate
        //what the water level should be.

        //We have that Volume = (Water Level)*(Cell Count)-(Total Elevation)
        //rearranging this gives us:
        water_level = (stdi.water_vol+total_elevation)/cells_affected.size();
      }


      //TODO: Use floating-point comparisons in these asserts.
      //Water level must be higher than (or equal to) the previous cell we looked at, but lower than (or equal to) the current cell
      std::cout<<"water level = "<<water_level<<" last topo "<<topo(cells_affected.back())<<" "<<bool(topo(cells_affected.back())<=water_level)<<std::endl;
      assert(cells_affected.size()==0 || topo(cells_affected.back())<=water_level);
      std::cout<<"water level = "<<water_level<<" my topo "<<topo(c.x,c.y)<<std::endl;
      assert(water_level<=topo(c.x,c.y));

      std::cerr<<"Adjusting wtd of depression...\n";
      std::cerr<<"\twater_level = "<<water_level<<std::endl;
      for(const auto c: cells_affected){
        std::cerr<<"Cell ("<<(c%topo.width())<<","<<(c/topo.width())<<") has elev="<<topo(c)<<", label="<<label(c)<<", wtd_old="<<wtd(c);
        assert(wtd(c)>=0);               //This should be true since we have been filling wtds as we go.
        assert(water_level>=topo(c));
        wtd(c) = water_level - topo(c);  //only change the wtd if it is an increase, here. We can't take water away from cells that already have it (ie reduce groundwater in saddle cells within a metadepression.)
        std::cerr<<", wtd_new="<<wtd(c)<<std::endl;
        assert(wtd(c)>=0);
      }

      //We've spread the water, so we're done        
      return;
      
    }  else {
      //We haven't found enough volume for the water yet.

      std::cout<<"we are in the else"<<std::endl;

      //During the adding of neighbours neighbours might get added that are
      //lower than we are and belong to a different depression (notably, this
      //happens at the edge of a flat abuting an ocean). These cells will then
      //be popped and could be processed inappropriately. To prevent this, we
      //skip them here.
      if(stdi.my_labels.count(label(c.x,c.y))==0){  //CHECK. This was preventing cells that flowed to the ocean from allowing my depression volume to update. Is this way ok? Is this even needed?
        std::cerr<<"Line "<<__LINE__<<std::endl;
        if(flood_q.empty())
          std::cout<<"that was the last cell"<<std::endl;
        continue;
      }

      //Okay, we're allow to add this cell's neighbours since this cell is part
      //of the metadepression.

      //Add this cell to those affected so that its volume is available for
      //filling.
      cells_affected.emplace_back(topo.xyToI(c.x,c.y));

      //Fill in cells' water tables as we go
      assert(wtd(c.x,c.y)<=0);
      stdi.water_vol += wtd(c.x,c.y);  //We use += because wtd is less than or equal to zero
      wtd(c.x,c.y)    = 0;             //Now we are sure that wtd is 0, since we've just filled it
      
      //Add the current cell's information to the running total
      total_elevation += topo(c.x,c.y);   //TODO: Should this be wtd? No. Since wtd is zero as of the lines just above.

      for(int n=0;n<neighbours;n++){
        std::cerr<<"in the for "<<n<<std::endl;
        const int nx = c.x + dx[n]; //TODO ModFloor(x+dx[n],topo.width()); //Get neighbour's x-coordinate using an offset and wrapping
        const int ny = c.y + dy[n];                     //Get neighbour's y-coordinate using an offset
        if(!topo.inGrid(nx,ny))   {                       //Is this cell in the grid?
          std::cerr<<"out of grid"<<std::endl;
          continue;                                     //Nope: out of bounds.
        }
     
        std::cerr<<"PQ inspecting ("<<nx<<","<<ny<<") which has label "<<label(nx,ny)<<std::endl;
        //Ocean cells may be found at the edge of a depression. They might get
        //added to this list even if there are other, lower, cells within the
        //depression which have not yet been explored. This happens when a flat
        //abutts an ocean. The side of the flat near the ocean will see the
        //ocean and try to add it. The ocean would then be called instead of
        //more cells within the depression. Therefore, we do not add ocean
        //cells.

        //We must use the ocean level rather than the ocean label, or we will
        //mistakenly miss adding higher cells which belong to the ocean's depression
        //e.g. an escarpment before the ocean. 

        if(!visited(nx,ny) && topo(nx,ny)>OCEAN_LEVEL){                  // add the neighbour only if it hasn't been added before 
          std::cerr<<"\tadding to the queue a value of "<<topo(nx,ny)<<" "<<" nx "<<nx<<" ny "<<ny<<" x "<<c.x<<" y "<<c.y<<std::endl;
          flood_q.emplace(nx,ny,topo(nx,ny));      //add all neighbours that haven't been added before to the queue. 
          visited(nx,ny) = true;
        }
      }
    }
  }

  //Since we're in this function we are supposed to be guaranteed to be able to
  //fill our depression, since we only enter this function once that is true.
  //Therefore, if we've reached this point, something has gone horribly wrong
  //somewhere. :-(

  std::cerr<<"PQ loop exited without filling a depression!"<<std::endl;
  
  std::cerr<<"Allowed labels = ";
  for(auto x:stdi.my_labels)
    std::cerr<<x<<" ";
  std::cerr<<std::endl;

  std::cerr<<"\tLabel of last cell       = "<<label(c.x,c.y)       <<std::endl;
  std::cerr<<"\tWater volume             = "<<stdi.water_vol       <<std::endl;
  std::cerr<<"\tCurrent volume           = "<<current_volume       <<std::endl;
  std::cerr<<"\tTotal elevation          = "<<total_elevation      <<std::endl;
  std::cerr<<"\tNumber of cells affected = "<<cells_affected.size()<<std::endl;  
  PrintDEM("Visited", visited);
  PrintDEM("Labels",  label  );
  throw std::runtime_error("PQ loop exited without filling a depression!");

}






  


template<class elev_t>
SubtreeDepressionInfo Find_filled(
  const int                               current_depression,    //Depression we are currently in
  const std::vector<Depression<elev_t> > &deps,                  //Depression hierarchy
  const rd::Array2D<float>               &topo,                  //Topographic data (used for determinining volumes as we're spreading stuff)
  const rd::Array2D<label_t>             &label,                 //Array indicating which leaf depressions each cell belongs to
  rd::Array2D<float>                     &wtd,                   //Water table depth
  std::string level =""                                          //TODO: For debugging
){
  //Stop when we reach one level below the leaves
  if(current_depression==NO_VALUE)
    return SubtreeDepressionInfo();

  std::cerr<<level<<"\033[93mInspecting depression "<<current_depression<<"\033[39m"<<std::endl;

  const auto& this_dep = deps.at(current_depression);

  //We start by visiting all of the ocean-linked depressions. They don't need to
  //pass us anything because their water has already been transferred to this
  //metadepression tree by Overflow(). Similar, it doesn't mater what their leaf
  //labels are since we will never spread water into them.
  for(const auto c: this_dep.ocean_linked)
    Find_filled(c, deps, topo, label, wtd, level+"\t");

  //At this point we've visited all of the ocean-linked depressions. Since all
  //depressions link to the ocean and the ocean has no children, this means we
  //have visited all the depressions and spread their water. Since we don't wish
  //to modify the ocean, we are done.
  if(current_depression==OCEAN)
    return SubtreeDepressionInfo();

  //We visit both of the children. We need to keep track of info from these
  //because we may spread water across them.
  SubtreeDepressionInfo left_info  = Find_filled(this_dep.lchild,deps,topo,label,wtd,level+"\t");
  SubtreeDepressionInfo right_info = Find_filled(this_dep.rchild,deps,topo,label,wtd,level+"\t");   

  SubtreeDepressionInfo combined;
  combined.my_labels.emplace(current_depression);
  combined.my_labels.merge(left_info.my_labels);
  combined.my_labels.merge(right_info.my_labels);

  //Add together marginal water volumes to produce a total water volume for the metadepression
  std::cerr<<"This water vol="<<this_dep.water_vol<<", LWV="<<left_info.water_vol<<", RWV="<<right_info.water_vol<<std::endl;
  combined.water_vol = this_dep.water_vol+left_info.water_vol+right_info.water_vol;

  combined.leaf_label = left_info.leaf_label;  //Choose left because right is not guaranteed to exist
  if(combined.leaf_label==NO_VALUE)            //If there's no label, then there was no child
    combined.leaf_label = current_depression;  //Therefore, this is a leaf depression

  combined.top_label = current_depression;



  std::cerr<<level<<"Total water volume is "<<combined.water_vol<<" and depression volume is "<<this_dep.dep_vol<<std::endl;

  //The water volume should never be greater than the depression volume because
  //otherwise we would have overflowed the water into the neighbouring
  //depression and moved the excess to the parent.
  if(combined.water_vol>this_dep.dep_vol){ //TODO: Make this an assert?
    throw std::runtime_error("water_vol>dep_vol");
  }



  //Since depressions store their marginal water volumes, if a parent depression
  //has 0 marginal water volume, then both of its children have sufficient
  //depression volume to store all of their water. However, if our parent is an
  //ocean-link then we are guaranteed to be able to fill now because excess
  //water will have been transferred into the parent and we don't want to pool
  //the parent's water with our own (it might be at the bottom of a cliff).

  if(deps.at(this_dep.parent).water_vol==0 || this_dep.ocean_parent){
    assert(combined.water_vol<=this_dep.dep_vol);

    //If both of a depression's children have already spread their water, we do not
    //want to attempt to do so again in an empty parent depression. 
    //We check to see if both children have finished spreading water. 

    Fill_Water(combined, deps, topo, label, wtd);

    //At this point there should be no more water all the way up the tree until
    //we pass through an ocean link, so we pass this up as a kind of null value.
    return SubtreeDepressionInfo();
  } else {
    return combined;
  }
}












  

     







  //TODO 4/5. Perform a depth-first post-order traversal of the depression
  //hierarchy (start with depressions for which `parent==NO_PARENT`. When you
  //reach the leaves if `water_vol>dep_vol` then try to overflow into the
  //neighbouring depression if its `water_vol<dep_vol`. To do so search
  //neighbour cells of `out_cell` for the lowest cell labeled `odep` and follow
  //that one's flow path until it terminates. Add excess water to that
  //depression's `water_vol`. After both child depressions have been visited
  //they will be finished trying to share their water and their excess water is
  //added to their parent's `water_vol`. Repeat the overflow attempt.

  //TODO 6: Adjust the hydrologic elevations. If a depression has `water_vol>0`
  //then all of its child depression are full. What remains is to use a
  //priority-queue and the Water Level Equation (see dephier.cpp) to find which
  //cells should be flooded. If a depression is a leaf depression (no children)
  //then start the PQ at the pit_cell. Otherwise, start at the pit cell of any
  //child depression since all children are flooded to at least the level of the
  //`out_elev` connecting the uppermost two children in the hierarchy. Cells
  //should be added to the queue with elevation equal to at least the `out_elev`
  //of the depression's children and, if `water_vol<dep_vol` should not exceed
  //the `out_elev` of the depression itself.





int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  const std::string in_name   = argv[1];
  const std::string out_name  = argv[2];
  const std::string out_graph = argv[3];

  rd::Array2D<float> topo = LoadData<float>(in_name,std::string("value"));   //Recharge (Percipitation minus Evapotranspiration)

  PrintDEM("Topography", topo);

  rd::Array2D<float>     wtd     (topo.width(), topo.height(), 1      ); //All cells have some water
  rd::Array2D<label_t>   label   (topo.width(), topo.height(), NO_DEP ); //No cells are part of a depression
  rd::Array2D<flowdir_t> flowdirs(topo.width(), topo.height(), NO_FLOW); //No cells flow anywhere

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++)
    if(topo(i)==OCEAN_LEVEL){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      label(i) = OCEAN;
      wtd  (i) = 0;
    }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = GetDepressionHierarchy<float,Topology::D8>(topo, label, flowdirs);

  PrintDEM("labels", label);

  std::cerr<<"\033[91m######################Depression Info\033[39m"<<std::endl;
  std::cerr<<std::setw(20)<<"Depression"<<std::setw(10)<<"Dep Vol"<<std::setw(10)<<"Water Vol"<<std::endl;
  for(unsigned int d=0;d<deps.size();d++)
    std::cerr<<std::setw(20)<<d<<std::setw(10)<<deps.at(d).dep_vol<<std::setw(10)<<deps.at(d).water_vol<<std::endl;
  std::cerr<<std::endl;

  //TODO: Remove. For viewing test cases.
  if(label.width()<1000){
  //  for(int y=0;y<label.height();y++){
      //for(int x=0;x<label.width();x++)
    //    std::cout<<std::setw(3)<<label(x,y)<<"checking ";
  //    std::cout<<std::endl;
//    }

    //GraphViz dot-style output for drawing depression hierarchy graphs.
    std::ofstream fgraph(out_graph);
    fgraph<<"digraph {\n";
    for(int i=0;i<(int)deps.size();i++){
      fgraph<<i<<" -> "<<deps[i].parent;
      if(deps[i].parent!=NO_VALUE && (deps[i].parent==OCEAN || !(deps[deps[i].parent].lchild==i || deps[deps[i].parent].rchild==i)))
        fgraph<<" [color=\"blue\"]";
      fgraph<<";\n";
    }
    fgraph<<"}\n";
  }

  SaveAsNetCDF(topo, out_name+"-topo.nc",       "value");
  SaveAsNetCDF(label,out_name+"-labels_raw.nc", "value");
  SaveAsNetCDF(label,out_name+"-labels_proc.nc","value");




  SurfaceWater(topo, wtd, label,deps,flowdirs);


  PrintDEM("wtd", wtd);





  Overflow(OCEAN, deps);

  std::cerr<<"\n";
  std::cerr<<std::setw(20)<<"Depression"<<std::setw(10)<<"Dep Vol"<<std::setw(10)<<"Water Vol"<<std::endl;
  for(unsigned int d=0;d<deps.size();d++)
    std::cerr<<std::setw(20)<<d<<std::setw(10)<<deps.at(d).dep_vol<<std::setw(10)<<deps.at(d).water_vol<<std::endl;
  std::cerr<<std::endl;

  PrintDEM("wtd", wtd, 9);

  for(auto &depression:deps){//(unsigned int d=0;d<deps.size();d++){
    std::cerr<<"Here's the list of all depressions with their parents and children: "
             <<depression.dep_label<<" "
             <<std::setw(3)<<depression.parent   <<" "
             <<std::setw(3)<<depression.lchild   <<" "
             <<std::setw(3)<<depression.rchild   <<" "
             <<depression.water_vol
             <<std::endl;
    std::cerr<<"\tOcean-linked = ";
    for(auto x: depression.ocean_linked)
      std::cerr<<x<<" ";
    std::cerr<<std::endl;
  }
  

  std::cerr<<"\n\n\033[91m#######################Finding Filled\033[39m"<<std::endl;





  Find_filled(OCEAN,deps,topo,label,wtd);                              //This should check everything that is an immediate child of the ocean, so we're supposed to hit all the depressions like this. 
       
  SaveAsNetCDF(wtd,out_name+"-wtd.nc","value");


  std::cerr<<"Finished"<<std::endl;

  if(topo.width()<1000){
    PrintDEM("topo",     topo    );
    PrintDEM("Flowdirs", flowdirs);
    PrintDEM("wtd",      wtd     );
    PrintDEM("labels",   label   );
  }

  return 0;
}
