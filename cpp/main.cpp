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
#include <richdem/flats/flats.hpp>
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
  for(int x=0;x<topo.width(); x++)
  for(int n=0;n<neighbours;n++){                           //Loop through neighbours
    const int nx = x+dx[n];                                //Identify coordinates of neighbour
    const int ny = y+dy[n];
    if(!topo.inGrid(nx,ny))
      continue;    
    if(flowdirs(nx,ny)==dinverse[n])  //CHECK WITH RICHARD           //(int)topo.xyToI(x,y))              //Does my neighbour flow into me?
      dependencies(x,y)++;                                 //Increment my dependencies
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
        std::cout<<"Giving water to "<<label(c)<<" "<<c<<" "<<deps[label(c)].water_vol<<std::endl;
   //   std::cout<<"the total water in this depression is "<<deps[label(c)].water_vol<<" and it is depression "<<label(c)<<std::endl;
        wtd(c) = 0; //Clean up as we go
      }
    } else {                               //not a pit cell
      //If we have water, pass it downstream.
      if(wtd(c)>0){ //Groundwater can go negative, so it's important to make sure that we are only passing positive water around
        std::cout<<"we have water coming from "<<c<<" and going to "<<n<<std::endl;
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

  //TODO: Does this really explore ALL the children?
  Overflow(this_dep.lchild, deps);

  //If this is false, but the previous was true, then we have found a depression that links to a depression that links to a ... that links to the ocean
  Overflow(this_dep.rchild, deps);

  //Catch depression that link to the ocean through this one. These are special
  //cases because we will never spread water across them and the current
  //depression: they only flow into the current depression

  for(const auto c: this_dep.ocean_linked)
    Overflow(c, deps);

  if(current_depression==OCEAN)
    return;

  std::cout<<"depression number "<<this_dep.dep_label<<" volume "<<this_dep.dep_vol<<" water "<<this_dep.water_vol<<std::endl;

  //Each depression will store the total volume of water it and all of its children contains.
  //That way the excess water a depression stores beyond that of its children is calculated as
  //Volume(Parent) - Volume(Children)

  //The depression is large enough to hold all of the water it contains, so we
  //stop trying to fill it. Since we are keeping track of total volumes we
  //propagate this water volume upward.
  if(this_dep.water_vol<=this_dep.dep_vol){
    deps.at(this_dep.parent).water_vol += this_dep.water_vol;
    return;
  }

  //We are overflowing the depression
  if(this_dep.odep == OCEAN){
    //If a depression overflows directly into an ocean then its odep is the
    //ocean and so is its parent.

    //The current depression's outlet is into the ocean. Since the ocean can
    //absorb an infinite amount of water without changing its water volume, we
    //simply set the amount of water contained in the current depression to be
    //equal to its depression volume. We throw away the excess water beyond this
    //since it does not affect the ocean.
    this_dep.water_vol = this_dep.dep_vol;
  } else {
   // std::cout<<"depression number "<<this_dep.dep_label<<" volume "<<this_dep.dep_vol<<" water "<<this_dep.water_vol<<std::endl;
    //The neighbouring depression is not the ocean; therefore, we should try to
    //overflow into it. To do so, we first calculate the excess water beyond
    //what this depression can hold (if this depression could hold all the
    //water, then we would have stopped above).
    float extra_water            = this_dep.water_vol - this_dep.dep_vol;

    //Now that we've figured out how much excess water there is, we mark this depression as being filled
    this_dep.water_vol = this_dep.dep_vol;

    auto &outlet_dep      = deps.at(this_dep.odep);              //Depression this one overflows into

    outlet_dep.water_vol += extra_water;                                   //Add excess water to the overflow depression
    if(outlet_dep.water_vol>outlet_dep.dep_vol){                           //Water in overflow depression exceeds its volume
      extra_water          = outlet_dep.water_vol - outlet_dep.dep_vol;    //Find out how much overflow depression overflowed by
      outlet_dep.water_vol = outlet_dep.dep_vol;                           //Mark overflow depression as being entirely filled
    }

    //Since if a depression links directly to an ocean its odep and parent
    //values are both OCEAN, we have overflowed into the ocean in the upper part
    //of this if-clause.
    assert(this_dep.parent!=OCEAN);

    //At this point extra_water is the amount of water that is left when both
    //this depression and its overflow depression are completely filled. We add
    //this extra water to this depression's (and therefore it's overflow
    //depression's) parent. We also add this depression's water volume since we
    //are tracking totals.
    deps.at(this_dep.parent).water_vol += this_dep.water_vol + extra_water;          //add any remaining water to the parent depression.
  
  }
  std::cout<<"and after: depression number "<<this_dep.dep_label<<" volume "<<this_dep.dep_vol<<" water "<<this_dep.water_vol<<std::endl;

}

 


//All overflowing depressions should by now have overflowed all the way down to the ocean. 
//We must now handle the actual water in the depressions and move it to wtd. 
//Completely full depressions and partially filled depressions must both update the wtd with their water. 



class SubtreeDepressionInfo {
 public:
  float water_vol = 0;
  float dep_vol   = 0;
  int   cells     = 0;
  int   bot_label = -1;
  int   top_label = -1;
  std::unordered_set<int> my_labels;
};




template<class elev_t>
void Fill_Water(
  SubtreeDepressionInfo                  &stdi,
  const std::vector<Depression<elev_t> > &deps,
  const rd::Array2D<float>               &topo,
  const rd::Array2D<label_t>             &label,
  rd::Array2D<float>                     &wtd
){
  //Nothing to do if we have no water
  if(stdi.water_vol==0)
    return;

  //changing tactics to start always from the leaves, then work your way up until you find something that isn't completely full. 
  std::cerr<<"\n\n\033[35m####################### Fill Water\033[39m"<<std::endl;
  rd::Array2D<bool> visited(topo.width(),topo.height(),false);

  const auto pit_cell    = deps.at(stdi.bot_label).pit_cell;
  double total_elevation = 0;
  std::cout<<"pit cell "<<pit_cell<<" "<<stdi.bot_label<<std::endl;

  GridCellZk_high_pq<elev_t> flood_q;                          

  std::cerr<<"Bottom label      = "<<stdi.bot_label<<std::endl;
  std::cerr<<"Depression volume = "<<stdi.dep_vol<<"\n";
  std::cerr<<"Water volume      = "<<stdi.water_vol<<std::endl;
  std::cerr<<"Allowed labels    = ";
  for(auto x:stdi.my_labels)
    std::cerr<<x<<" ";
  std::cerr<<std::endl;

  assert(pit_cell>=0);

  flood_q.emplace(
    pit_cell % topo.width(),
    pit_cell / topo.width(),
    topo(pit_cell)
  );                    //create a new priority queue starting with the pit cell of the depression

  visited(pit_cell) = true;//label(pit_cell);         // show that we have already added this cell to those that have water. We need a better way to do this. 

  double current_volume;
  GridCellZk_high<elev_t> c(0,0,0,0);

  std::vector<int> cells_affected;

  std::cout<<"about to start the flood_q while loop "<<label(1,1)<<std::endl;
  while(!flood_q.empty()){
    c = flood_q.top(); //TODO local var
    flood_q.pop();

    std::cout<<"width "  <<topo.width()<<std::endl;
    std::cout<<"x and y "<<c.x<< " "<<c.y<<std::endl;
    std::cout<<"current cell elevation "<<topo(c.x,c.y)<<std::endl;

    //We have two volumes to keep track of. One is the above ground volume formed by the virtual waterline we are raising,
    //the other is the volume of water that the water table can absorb. We track `total_wtd` as a summed variable since
    //we will always fill it, if possible. We keep track of above ground volume as a multiplied variable since we will
    //adjust the final water level only above the ground.

    //TODO: Note that the current cell's above ground volume and wtd do not contribute at all. This a choice that Kerry and Richard discussed
    //It is as though there is a virtual water line coincident with the edge of the current cell. No water infiltrates into this cell
    //or is stored above it - only cells previously visited are considered when doing volume calculations.
    current_volume = (cells_affected.size()*topo(c.x,c.y) - total_elevation);       //get the current volume of this part of the depression //TODO: local var

    std::cout<<"volume "<<current_volume<<" "<<cells_affected.size()<<" "<<topo(c.x,c.y)<<" "<<total_elevation<<std::endl;

    if(stdi.water_vol<0 && stdi.water_vol>-1e-6) //TODO: Sprinkle this liberally everywhere
      stdi.water_vol = 0;
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

    if(stdi.water_vol<=current_volume-wtd(c.x,c.y)){         //if this volume will accommodate all of the water, stop adding cells
      const auto my_elev = topo(c.x,c.y);

      std::cerr<<"Attempting to fill depression..."<<std::endl;
      std::cerr<<"\tLabel of last cell       = "<<label(c.x,c.y)       <<std::endl;
      std::cerr<<"\tWater volume             = "<<stdi.water_vol       <<std::endl;
      std::cerr<<"\tDepression volume        = "<<stdi.dep_vol         <<std::endl;
      std::cout<<"\tDepression number        = "<<stdi.bot_label       <<std::endl;
      std::cerr<<"\tCurrent volume           = "<<current_volume       <<std::endl;
      std::cerr<<"\tTotal elevation          = "<<total_elevation      <<std::endl;
      std::cerr<<"\tCurrent elevation        = "<<my_elev              <<std::endl;
      std::cerr<<"\tNumber of cells affected = "<<cells_affected.size()<<std::endl;  

      //Above we only required that the volume of water be less than the current
      //volume plus the available infiltration space of the current cell. If the
      //current volume is not enough to hold all the water, that we means we
      //have to stash some of it in the current cell's water table.
      if(current_volume<stdi.water_vol){ //TODO: Check stdi.my_labels.count(label(c.x,c.y))==0 ?
        //Fill in as much of this cell's water table as we can
        const double fill_amount = stdi.water_vol - current_volume;
        assert(fill_amount>=0);
        wtd(c.x,c.y)   += fill_amount;
        stdi.water_vol -= fill_amount;
      }

      //At this point we know that the surface of the water should be at or below the elevation of the cell
      const double water_level = (stdi.water_vol+total_elevation)/cells_affected.size();

      //Water level must be higher than (or equal to) the previous cell we looked at, but lower than (or equal to) the current cell
      assert(cells_affected.size()==0 || topo(cells_affected.back())<=water_level);
      assert(water_level<=topo(c.x,c.y));

      std::cerr<<"Adjusting wtd of depression...\n";
      std::cerr<<"\twater_level = "<<water_level<<std::endl;
      for(const auto c: cells_affected){
        std::cerr<<"Cell ("<<(c%topo.width())<<","<<(c/topo.width())<<") has elev="<<topo(c)<<", label="<<label(c)<<", wtd_old="<<wtd(c);
        wtd(c) = water_level - topo(c);
        std::cerr<<", wtd_new="<<wtd(c)<<std::endl;
        assert(wtd(c)>=0);
      }
        
      return;
      
    } else {
      //During the adding of neighbours neighbours might get added are are lower than we are and belong to a different
      //depression (notably, this happens at the edge of a flat abuting an ocean). These cells will then be popped
      //and could be processed inappropriately. To prevent this, we skip them here.
      if(stdi.my_labels.count(label(c.x,c.y))==0)
        continue;

      //We haven't found enough volume for the water yet. Add this cell's neighbours and move on to the next
      //highest cell, which will add the volume of the cell we've just visited

      //Add the cell we've just visited, so that its volume is available in the next iteration
      cells_affected.emplace_back(topo.xyToI(c.x,c.y));   //list all of the cells in which the wtd will be changing

      //Fill in cells' water tables as we go
      stdi.water_vol += wtd(c.x,c.y);   
      wtd(c.x,c.y)    = 0;
      
      //Add the current cell's information to the running total
      total_elevation += topo(c.x,c.y);   //TODO: Should this be wtd?

      //TODO: Use labels to prevent visiting places we shouldn't
      for(int n=0;n<neighbours;n++){
        const int nx = c.x + dx[n]; //TODO ModFloor(x+dx[n],topo.width()); //Get neighbour's x-coordinate using an offset and wrapping
        const int ny = c.y + dy[n];                     //Get neighbour's y-coordinate using an offset
        if(!topo.inGrid(nx,ny))                          //Is this cell in the grid?
          continue;                                     //Nope: out of bounds.
     
        std::cout<<"PQ inspecting ("<<nx<<","<<ny<<") which has label "<<label(nx,ny)<<std::endl;
        //Ocean cells may be found at the edge of a depression. They might get
        //added to this list even if there are other, lower, cells within the
        //depression which have not yet been explored. This happens when a flat
        //abutts an ocean. The side of the flat near the ocean will see the
        //ocean and try to add it. The ocean would then be called instead of
        //more cells within the depression. Therefore, we do not add ocean
        //cells.
        if(!visited(nx,ny) && topo(nx,ny)>OCEAN_LEVEL){                  // add the neighbour only if it hasn't been added before 
          std::cout<<"\tadding to the queue a value of "<<topo(nx,ny)<<" "<<" nx "<<nx<<" ny "<<ny<<" x "<<c.x<<" y "<<c.y<<std::endl;
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






  


//TODO: Use a hashset to keep track of which depression labels are valid for filling
template<class elev_t>
SubtreeDepressionInfo Find_filled(
  const int                               current_depression,
  const std::vector<Depression<elev_t> > &deps, 
  const rd::Array2D<float>               &topo,  
  const rd::Array2D<label_t>             &label, 
  rd::Array2D<float>                     &wtd,
  std::string level ="" //TODO
){
  //Stop when we reach one level below the leaves

  if(current_depression==NO_VALUE)
    return SubtreeDepressionInfo();

  std::cerr<<level<<"\033[93mInspecting depression "<<current_depression<<"\033[39m"<<std::endl;

  const auto& this_dep = deps.at(current_depression);

  SubtreeDepressionInfo left_info  = Find_filled(this_dep.lchild,deps,topo,label,wtd,level+"\t");                              //This should check everything that is an immediate child of the ocean, so we're supposed to hit all the depressions like this. 
  SubtreeDepressionInfo right_info = Find_filled(this_dep.rchild,deps,topo,label,wtd,level+"\t");   

  for(const auto c: this_dep.ocean_linked)
    SubtreeDepressionInfo temp =  Find_filled(c, deps, topo, label, wtd, level+"\t"); //TODO: need to add to combined

  if(current_depression==OCEAN)
    return SubtreeDepressionInfo();              //better place to put this? CHECK

  std::cerr<<level<<"Water volume is "<<this_dep.water_vol<<" and depression volume is "<<this_dep.dep_vol<<std::endl;

  //The water volume should never be greater than the depression volume because
  //otherwise we would have overflowed the water into the neighbouring
  //depression and moved the excess to the parent.
  if(this_dep.water_vol>this_dep.dep_vol){ //TODO: Make this an assert?
    throw std::runtime_error("water_vol>dep_vol");
  }

  SubtreeDepressionInfo combined;
  combined.water_vol = this_dep.water_vol;//left_info.water_vol + right_info.water_vol + this_dep.water_vol; //Water volumes are now total, not marginal
  std::cerr<<level<<"Getting combined water volumes "<<left_info.water_vol<<" "<<right_info.water_vol<<" "<<this_dep.water_vol<<std::endl;
  
  combined.dep_vol   = this_dep.dep_vol;                                                //Depression volumes are total
  std::cerr<<level<<"Getting combined depression volumes "<<left_info.dep_vol<<" "<<right_info.dep_vol<<" "<<this_dep.dep_vol<<std::endl;
  
  combined.cells     = left_info.cells + right_info.cells;

  combined.my_labels.emplace(current_depression);
  combined.my_labels.merge(left_info.my_labels);
  combined.my_labels.merge(right_info.my_labels);

  if(left_info.bot_label==NO_VALUE)
    combined.bot_label = current_depression;
  else
    combined.bot_label = left_info.bot_label; //Choose left because right is not guaranteed to exist

  combined.top_label = current_depression;

  //Depressions that overflow have water_vol=dep_vol, so we don't want to spread
  //that water around yet. (There's a special case in which a depression has
  //*exactly* water_vol=dep_vol, but this is taken care of by its parent)
  //However, once we've reached a point where our parent is ocean-linking, we've
  //reached the farthest point that water can be passed up and need to spread it
  //around.

  //Why is this so? Because if we go up farther we will either end up spreading
  //water across the ocean, or across an ocean-linking depression. In general,
  //ocean-linking depressions are BELOW us (like at the bottom of a tall cliff),
  //so trying to spread across this depression and its parent doesn't make
  //sense.

  //We want to avoid spreading water around more than once. Therefore, we want
  //to find the largest grouping of depressions which will share water. Thus,
  //when we see `combined.water_vol==combined.dep_vol`, this means we *could*
  //spread water around, but that would waste time since we might spread the
  //metadepression's water across the same cells. Therefore, we skip cases where
  //`water_vol==dep_vol` and only trigger on cases where `water_vol<dep_vol`.

  //However, if the parent of this depression is the ocean-linking, then we
  //don't want to try spreading water across the whole ocean (or something at
  //the bottom of a cliff) as well as the metadepression. Therefore, if this
  //depression's parent is ocean-linking, we stop and spread the water, even if
  //`water_vol==dep_vol`

  //There is another case in which we want to immediately spread water. Imagine
  //the situation below:
  //      ___ _ _ _  _C_       ___
  //         \  A   /   \  B  /
  //          \    /     \_ _/
  //           \__/       \_/
  //We are in depression A and its water volume exactly equals its depression
  //volume, so we should spread. But we can't capture this using
  //`water_vol<=dep_vol` for reasons explained above. However, we need to spread
  //now because if we delegate spreading to C, then C might place all of A's
  //water into B, which would not be appropriate. Therefore, if
  //`water_vol==dep_vol` and B's `water_vol<dep_vol`, then B has not overflowed
  //into C so we should be spreading A and B separately, rather than adding them
  //both to C and letting C spread across both of them.

  const auto odep = deps.at(deps.at(combined.top_label).odep);

  if(combined.water_vol<combined.dep_vol || this_dep.ocean_parent || (combined.water_vol==combined.dep_vol && odep.water_vol!=odep.dep_vol)){
    //If this depression's parent is the ocean, we want to spread the water now.
    //If we wait, we'll try to spread it across the whole ocean. However, we
    //have to guarantee that there is enough volume in the depression to contain
    //the water.
    if(this_dep.parent==OCEAN && combined.water_vol>combined.dep_vol){
      std::cerr<<"this is dep "<<this_dep.dep_label<<" vol "<<this_dep.dep_vol<<" "<<" water "<<this_dep.water_vol<<std::endl;
      std::cerr<<"combined vol "<<combined.dep_vol<<" "<<" water "<<combined.water_vol<<std::endl;
      throw std::logic_error("Water volume must be <= depression volume if depression's parent is the ocean.");
    }

    //If both of a depression's children have already spread their water, we do not
    //want to attempt to do so again in an empty parent depression. 
    //We check to see if both children have finished spreading water. 

    if(    (this_dep.lchild!=NO_VALUE && left_info.water_vol>0)                      //TODO: This would be easier to deal with if we switch back to marginal water volume
        || (this_dep.rchild!=NO_VALUE && right_info.water_vol>0) 
        || (this_dep.lchild == NO_VALUE && this_dep.rchild == NO_VALUE)
    ){
      Fill_Water(combined, deps, topo, label, wtd);
    }

    SubtreeDepressionInfo temp;
    return temp;
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

  // richdem::ResolveFlatsEpsilon(topo); //TODO: Shouldn't need this.

  rd::Array2D<float> wtd (topo.width(),topo.height(),1);

 
  //Initialize labels to indicate that none of the cells are part of depressions
  rd::Array2D<label_t> label   (topo.width(),topo.height(),NO_DEP);

  //Initialize flow directions to indicate that none of the cells flow anywhere
  flowdirs.resize(topo.width(),topo.height(),NO_FLOW);

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++)
    if(topo(i)==OCEAN_LEVEL)
      label(i) = OCEAN;

  //Label all the depressions and get the hierarchy connecting them.

  auto deps = GetDepressionHierarchy<float,Topology::D8>(topo, label, flowdirs);

 for(auto &depression: deps){
    std::cout<<" we are on depression "<<depression.dep_label<<std::endl;
    std::cout<<"all the pit cells are "<<depression.pit_cell<<std::endl;
 }



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

  SaveAsNetCDF(topo,out_name+"-topo.nc","value");
  SaveAsNetCDF(label,out_name+"-labels_raw.nc","value");

  // LastLayer(label, topo, deps);

  SaveAsNetCDF(label,out_name+"-labels_proc.nc","value");



  //Set all ocean cells to 0 wtd
  for(unsigned int i=0;i<wtd.size();i++){
    wtd(i)=1;
    if(topo(i)<=OCEAN_LEVEL)
      wtd(i)=0;
  }

  //for(auto &depression: deps){
//std::cout<<depression.lchild<<"label "<<depression.dep_label<<std::endl;
//}

  SurfaceWater(topo, wtd, label,deps,flowdirs);

  PrintDEM("labels", label);
  PrintDEM("wtd", wtd);

  for(int y=0;y<topo.height();y++){                       //add all of the cells
    for(int x=0;x<topo.width();x++){
      std::cout<<"x "<<x<<" y "<<y<<" flow "<<int(flowdirs(x,y))<<std::endl;
    }
  }


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



 // Fill_Water(deps,topo,label,wtd);
       
  SaveAsNetCDF(wtd,out_name+"-wtd.nc","value");


  std::cout<<"LIne "<<__LINE__<<std::endl;

 // std::ofstream fout("/z/out.dat");

 // fout.write(reinterpret_cast<const char*>(wtd.data), wtd.size()*sizeof(float));

  if(topo.width()<1000){
    PrintDEM("topo", topo, 9);
    PrintDEM("Flowdirs", flowdirs, 9);
    PrintDEM("wtd", wtd, 9);
  }

  return 0;
}

