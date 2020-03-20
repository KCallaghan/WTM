#ifndef _df_flow_hpp_
#define _df_flow_hpp_

#include "dephier.hpp"
#include "DisjointDenseIntSet.hpp"
#include "../common/netcdf.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/depressions/depressions.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace richdem::dephier {

const int *const dx         = d8x;
const int *const dy         = d8y;
const int *const dinverse   = d8_inverse;
const int        neighbours = 8;

const double FP_ERROR = 1e-4;
const double infiltration_FSM = 0.5;
//const double evaporation_FSM = 0.5; //TODO - obviously this should be based on actual surface evap rates. This is just for testing purposes quickly.


///////////////////////////////////
//Function prototypes
///////////////////////////////////

template<class elev_t, class wtd_t>
void FillSpillMerge(
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<dh_label_t> &label,
  const rd::Array2D<dh_label_t> &final_label,
  const rd::Array2D<flowdir_t>  &flowdirs,
  DepressionHierarchy<elev_t>   &deps,
  rd::Array2D<wtd_t>            &wtd,
  ArrayPack                     &arp
);

template<class elev_t, class wtd_t>
static void MoveWaterIntoPits(
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<int>        &label,
  const rd::Array2D<dh_label_t> &final_label,
  const rd::Array2D<flowdir_t>  &flowdirs,
  DepressionHierarchy<elev_t>   &deps,
  rd::Array2D<wtd_t>            &wtd,
  ArrayPack                     &arp
);


template<class elev_t, class wtd_t>
static void CalculateWtdVol(
  rd::Array2D<wtd_t>            &wtd,
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<dh_label_t> &final_label,
  DepressionHierarchy<elev_t>   &deps
);


template<class elev_t,class wtd_t>
static void MoveWaterInDepHier(
  int                                        current_depression,
  DepressionHierarchy<elev_t>                &deps,
  const rd::Array2D<elev_t>                  &topo,
  const rd::Array2D<int>                     &label,
  const rd::Array2D<int>                     &final_label,
  const rd::Array2D<flowdir_t>               &flowdirs,
  rd::Array2D<wtd_t>                         &wtd,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  ArrayPack                                  &arp
);


template<class elev_t,class wtd_t>
static void MoveWaterInOverflow(
  float                          extra_water,
  const dh_label_t               current_dep,
  const dh_label_t               previous_dep,
  const rd::Array2D<elev_t>      &topo,
  rd::Array2D<wtd_t>             &wtd,
  const rd::Array2D<flowdir_t>   &flowdirs,
  const rd::Array2D<int>         &final_label,
  const rd::Array2D<int>         &label,
  DepressionHierarchy<elev_t>    &deps,
  ArrayPack                      &arp
  );

template<class elev_t, class wtd_t>
static dh_label_t OverflowInto(
  const dh_label_t                           root,
  const dh_label_t                           previous_dep,
  const dh_label_t                           stop_node,
  const rd::Array2D<elev_t>                  &topo,
  const rd::Array2D<int>                     &label,
  const rd::Array2D<int>                     &final_label,
  const rd::Array2D<flowdir_t>               &flowdirs,
  rd::Array2D<wtd_t>                         &wtd,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  double                                     extra_water,
  ArrayPack                                  &arp
);


class SubtreeDepressionInfo;
template<class elev_t>
static SubtreeDepressionInfo FindDepressionsToFill(
  const int                         current_depression, 
  const DepressionHierarchy<elev_t> &deps,               
  const rd::Array2D<float>          &topo,               
  const rd::Array2D<dh_label_t>     &label,              
  rd::Array2D<float>                &wtd,           
  ArrayPack                         &arp     
);

template<class elev_t, class wtd_t>
static void FillDepressions(
  SubtreeDepressionInfo             &stdi,  
  double                            water_vol, 
  const DepressionHierarchy<elev_t> &deps,      
  const rd::Array2D<float>          &topo,      
  const rd::Array2D<dh_label_t>     &label,     
  rd::Array2D<wtd_t>                &wtd,   
  ArrayPack                         &arp     
);



///////////////////////////////////
//Implementations
///////////////////////////////////



///This function routes surface water from into pit cells and then distributes
///it so that it fills the bottoms of depressions, taking account of overflows.
///
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression 
///                each cell belongs to.
///@param flowdirs Flowdirs generated by GetDepressionHierarchy
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param wtd      Water table depth. Values of 0 indicate saturation. 
///                Negative values indicate additional water can be added to the
///                cell. Positive values indicate standing surface water.
///
///Note that from GetDepressionHierarchy we already know the number of cells
///within each depression as well as its volume.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water contained in each depression. Modifies `wtd` to indicate how
///        saturated a cell is or how much standing surface water it has.
template<class elev_t, class wtd_t>
void FillSpillMerge(
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<dh_label_t> &label,
  const rd::Array2D<dh_label_t> &final_label,
  const rd::Array2D<flowdir_t>  &flowdirs,
  DepressionHierarchy<elev_t>   &deps,
  rd::Array2D<wtd_t>            &wtd,
  ArrayPack                     &arp
){
  rd::Timer timer_overall;
  timer_overall.start();
  
  //We move standing water downhill to the pit cells of each depression
  MoveWaterIntoPits(topo, label,final_label, flowdirs, deps, wtd,arp);

  { 
    //Scope to limit `timer_overflow` and `jump_table`. Also ensures
    //`jump_table` frees its memory
    rd::Timer timer_overflow;
    timer_overflow.start();
    std::unordered_map<dh_label_t, dh_label_t> jump_table;
 
    //calculate the wtd_vol of depressions, in order to be able to know which need to overflow and which can accommodate more water:
    //the reason for only calculating it here is that it's better to do this after MoveWaterIntoPits, when a lot of the infiltration
    //that would occur has already happened, so we didn't have to constantly update wtd_vol during MoveWaterIntoPits. 
    CalculateWtdVol(wtd,topo,final_label,deps,arp);
 
    //Now that the water is in the pit cells, we move it around so that
    //depressions which contain too much water overflow into depressions that
    //have less water. If enough overflow happens, then the water is ultimately
    //routed to the ocean.
    MoveWaterInDepHier(OCEAN, deps, topo,label,final_label,flowdirs,wtd,jump_table,arp);

    std::cerr<<"t FlowInDepressionHierarchy: Overflow time = "<<timer_overflow.stop()<<std::endl;
  }

  //Sanity checks
  for(int d=1;d<(int)deps.size();d++){
    const auto &dep = deps.at(d);
 
    assert(dep.water_vol==0 || dep.water_vol<=dep.wtd_vol);
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.lchild!=NO_VALUE && deps.at(dep.lchild).water_vol - dep.water_vol < FP_ERROR));
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.rchild!=NO_VALUE && deps.at(dep.rchild).water_vol - dep.water_vol < FP_ERROR));
    if(dep.lchild != NO_VALUE && deps.at(dep.lchild).water_vol > dep.water_vol)
      deps.at(dep.lchild).water_vol = dep.water_vol;
    if(dep.rchild != NO_VALUE && deps.at(dep.rchild).water_vol > dep.water_vol)
      deps.at(dep.rchild).water_vol = dep.water_vol;

  }

  std::cerr<<"p Finding filled..."<<std::endl;
  rd::Timer timer_filled;
  timer_filled.start();
  //We start at the ocean, crawl to the bottom of the depression hierarchy and
  //determine which depressions or metadepressions contain standing water. We
  //then modify `wtd` in order to distribute this water across the cells of the
  //depression which will lie below its surface.
  FindDepressionsToFill(OCEAN,deps,topo,label,wtd,arp);                              
  std::cerr<<"t FlowInDepressionHierarchy: Fill time = "<<timer_filled.stop()<<" s"<<std::endl;

  std::cerr<<"t FlowInDepressionHierarchy = "<<timer_overall.stop()<<" s"<<std::endl;
}



///This function works much like a standard flow accumulation algorithm except
///that as the water moves downhill it contributes to saturating the water table
///and, when it reaches the pit cell of a depression, all the remaining water is
///moved into the DepressionHierarchy data structure for rapid flood-spill-merge
///calculations.
///
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression 
///                each cell belongs to.
///@param flowdirs Flowdirs generated by GetDepressionHierarchy
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param wtd      Water table depth. Values of 0 indicate saturation. 
///                Negative values indicate additional water can be added to the
///                cell. Positive values indicate standing surface water.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water contained in each LEAF depression. Modifies `wtd` to indicate 
///        how saturated a cell is. All values in wtd will be <=0 following this
///        operation.
template<class elev_t, class wtd_t>
static void MoveWaterIntoPits(
  const rd::Array2D<elev_t>    &topo,
  const rd::Array2D<int>       &label,
  const rd::Array2D<int>       &final_label,
  const rd::Array2D<flowdir_t> &flowdirs,
  DepressionHierarchy<elev_t>  &deps,
  rd::Array2D<wtd_t>           &wtd,
  ArrayPack                    &arp
){
  rd::Timer timer;
  rd::ProgressBar progress;
  timer.start();

  //Our first step is to move all of the water downstream into pit cells. To do
  //so, we use the steepest-descent flow directions provided by the depression
  //hierarchy code

  std::cerr<<"p Moving surface water downstream..."<<std::endl;


  #pragma omp parallel for 
  for(unsigned int i=0;i<topo.size();i++){   
    if(wtd(i)>0){
      arp.surface_water(i) += wtd(i);    //Move any surface water into a separate array. This is to allow us to have infiltration into wtd, that does not necessarily fill all the way to the surface, and still have water continue moving downslope. 
      wtd(i) = 0;                        //And we reset any cells that contained surface water to 0. 
    }
  }
 

  //Calculate how many upstream cells flow into each cell
  rd::Array2D<char>  dependencies(topo.width(),topo.height(),0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width(); x++)
  for(int n=1;n<=neighbours;n++){     //Loop through neighbours
    const int nx = x+dx[n];           //Identify coordinates of neighbour
    const int ny = y+dy[n];
    if(!topo.inGrid(nx,ny))
      continue;    
    if(flowdirs(nx,ny)==dinverse[n])  //Does my neighbour flow into me?
      dependencies(x,y)++;            //Increment my dependencies
  }

  
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

  for(int d=0;d<(int)deps.size();d++){   //reset all of the water volumes to 0 for each time we iterate through the surface water. 
    auto &dep = deps.at(d);                
    dep.water_vol = 0;
    dep.wtd_vol = 0;
  }

  
  //Starting with the peaks, pass flow downstream
  progress.start(topo.size());
  while(!q.empty()){
    ++progress;

    const auto c = q.front();          //Copy focal cell from queue
    q.pop();                           //Clear focal cell from queue

    //Coordinates of downstream neighbour, if any
    const auto ndir = flowdirs(c); 
    int x,y,nx,ny;
    topo.iToxy(c,x,y);

    int n = NO_FLOW;
    if(ndir!=NO_FLOW){  //TODO: Fix this monkey patching
      nx = x+dx[ndir];
      ny = y+dy[ndir];
      n            = topo.xyToI(nx,ny);
      assert(n>=0);
    }


    if(label(c) == OCEAN){    //If downstream neighbour is the ocean, we drop our water into it and the ocean is unaffected. 
      arp.surface_water(c) = 0;
    }
  

    if(n == NO_FLOW){    //if this is a pit cell, move the water to the appropriate depression's water_vol.    
      if(arp.surface_water(c)>0){
        deps[label(c)].water_vol += arp.surface_water(c)*arp.cell_area[y];   //surface_water is a depth of water and water_vol is a volume, so multiply by the area of the pit cell to convert. 
        assert(deps[label(c)].water_vol >= -FP_ERROR);
        if(deps[label(c)].water_vol < 0)
          deps[label(c)].water_vol = 0.0;
        arp.surface_water(c) = 0; //Clean up as we go
      }
    }else {                               //not a pit cell

//let some evaporation happen
      if(arp.surface_water(c)>0){
        arp.surface_water(c) = arp.surface_water(c)*(1-surface_evap(c));
      }

      //If we still have water, pass it downstream.
      if(arp.surface_water(c)>0){ 
        //We use cell areas because the depth will change if we are moving to a cell on a different latitude. 
        arp.surface_water(n) += arp.surface_water(c)*arp.cell_area[y]/arp.cell_area[ny];  //Add water to downstream neighbour. This might result in filling up the groundwater table, which could have been negative
        arp.surface_water(c)  = 0;       //Clean up as we go

        if(wtd(n)<0){  //This is where we allow water to infiltrate, but not all of it. Only some proportion, dictated by the value of infiltration_FSM. 
          float infiltration_amount = arp.surface_water(n)*infiltration_FSM;
          if(arp.surface_water(n)<=0.001)
            infiltration_amount = arp.surface_water(n);     //If there is a mm or less of surface water, just let it all infiltrate. Else we end up with super small bits of water still moving around and it's just silly later on. 
          if((wtd(n) + infiltration_amount) <= 0){
            wtd(n) += infiltration_amount;
            arp.surface_water(n) -= infiltration_amount; 
          }
          else{
            arp.surface_water(n) += wtd(n); //plus because wtd is negative. 
            wtd(n) = 0;
          }
        }

      }



      //Decrement the neighbour's dependencies. If there are no more dependencies,
      //we can process the neighbour.
      if(--dependencies(n)==0){                            
        assert(dependencies(n)>=0);
        q.emplace(n);                   //Add neighbour to the queue
      }
    }
  }

  progress.stop();

  std::cerr<<"t FlowInDepressionHierarchy: Surface water = "<<timer.stop()<<" s"<<std::endl;
}




///Here, we calculate the wtd_vol of all of the depressions. 
///The reason for doing it here rather than in dephier is that I avoid 
///having to make many adjustments to the wtd_vols during MoveWaterIntoPits.
///Wtd_vol represents the dep_vol plus any additional volume available to 
///store water as groundwater. 
///
///
///@param wtd      Water table depth. Values of 0 indicate saturation. 
///                Negative values indicate additional water can be added to the
///                cell. Positive values indicate standing surface water.
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression 
///                each cell belongs to.
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water storage space available, including groundwater storage.
template<class elev_t, class wtd_t>
static void CalculateWtdVol(
  rd::Array2D<wtd_t>            &wtd,
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<dh_label_t> &final_label,
  DepressionHierarchy<elev_t>   &deps,
  ArrayPack                     &arp
){

  for(int d=0;d<(int)deps.size();d++){
    auto &dep = deps.at(d);

    if(dep.dep_label==OCEAN)
      continue;
    dep.wtd_only = 0;    //reset all of the wtd heights so that they don't just get added from previous iterations. 
    dep.wtd_vol = dep.dep_vol;  //We start off with the depression volume and then just have to add the additional below-ground volume to this. 
    }

  //We need to create yet another measure of volume which I 
  //for now will call wtd_vol. This will be the dep_vol plus the additional
  //volume allowed in the depression through storage as groundwater. 
  //This is important because a depression may actually be able to store more 
  //water than in its dep_vol. We may otherwise be overflowing a depression when it
  //is not actually supposed to overflow! 
  //When we do overflow, we will also have to keep track of changes to the wtd
  //in the overflow depression, and associated changes in wtd_vol. When a depression
  //is completely saturated in groundwater, we will have wtd_vol == dep_vol.

  for(int y=0;y<wtd.height();y++)
  for(int x=0;x<wtd.width();x++){  //cycle through the domain and add up all of the below-ground water storage space available
    auto clabel        = final_label(x,y);    

    if(clabel==OCEAN)
      continue;

    assert(wtd(x,y) <= FP_ERROR);  //because all wtds get set to 0 when doing movewaterintopits - surface water is now gathered in the pits and has not yet been distributed.
    if(wtd(x,y) > 0)
      wtd(x,y) = 0.0;

    deps[clabel].wtd_only -= wtd(x,y)*arp.cell_area[y];  //- because wtd is negative. This records only the below-ground volume available
    deps[clabel].wtd_vol  -= wtd(x,y)*arp.cell_area[y];  //and this records the total volume - above and below ground - available to store water. 
    }


  for(int d=0;d<(int)deps.size();d++){
    auto &dep = deps.at(d);
    if(dep.dep_label==OCEAN)
      continue;
    if(dep.lchild!=NO_VALUE){
  
      dep.wtd_vol      += deps.at(dep.lchild).wtd_only;         //store the total wtd_vols with all of your children included. 
      dep.wtd_vol      += deps.at(dep.rchild).wtd_only;

      dep.wtd_only      += deps.at(dep.lchild).wtd_only;  
      dep.wtd_only      += deps.at(dep.rchild).wtd_only;  

    }
  
    assert(dep.wtd_vol >= -FP_ERROR);
    if(dep.wtd_vol < 0)
      dep.wtd_vol = 0.0;
    assert(dep.wtd_vol - dep.dep_vol >= -FP_ERROR);
    if(dep.wtd_vol < dep.dep_vol)
      dep.wtd_vol = dep.dep_vol;
  }


for(int d=0;d<(int)deps.size();d++){  //Just checking that nothing went horribly wrong. 
    auto &dep = deps.at(d);
    if(dep.dep_label==OCEAN)
      continue;
    if(dep.lchild!=NO_VALUE){
      assert(dep.dep_vol>= deps.at(dep.lchild).dep_vol);
      assert(dep.dep_vol>= deps.at(dep.lchild).dep_vol);
      assert(dep.wtd_vol>= deps.at(dep.lchild).wtd_vol);
      assert(dep.wtd_vol>= deps.at(dep.rchild).wtd_vol);
    }
  }
}



///At this point all the values of the water table `wtd`, which is not used by
///this function, are <=0. The excess water, which will eventually become
///standing surface water, is stored in the DepressionHierarchy itself in the
///leaf depressions.
///
///In this function we will perform a depth-first post-order traversal of the
///depression hierarchy starting with the OCEAN. When we reach the leaf
///depressions we check if `water_vol>dep_vol`. If so, we try to overflow into
///the geographically proximal leaf depression indicaed by our outlet. If there
///is not sufficient room in the depression linked to by our outlet to hold all
///the water, then the excess is passed to the depression's parent.
///
///Thus, by the time we exit this function, water will have been redistributed
///from leaf depressions as far up the depression hierarchy as necessary to
///ensure that there is sufficient volume to hold it. Any excess water is
///redistributed to the OCEAN, which is unaffected.
///
///The modified depression hierarchy's values for depression and water volume
///indicate total volumes: the volume in just the metadepression as well as its
///children.
///
///@param current_depression  We're doing a tree traversal; this is the id of 
///                           the depression we're currently considering.
///@param deps                The DepressionHierarchy generated by 
///                           GetDepressionHierarchy
///@param jump_table          A data structure that persists throughout the
///                           traversal and is used to skip from leaf nodes to
///                           the highest known meta-depression which still has
///                           unfilled volume. Ensures the traversal happens in
///                           O(N) time.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water in each depression. This information can be used to add
///        standing surface water to the cells within a depression.
template<class elev_t,class wtd_t>
static void MoveWaterInDepHier(
  int                                        current_depression,
  DepressionHierarchy<elev_t>                &deps,
  const rd::Array2D<elev_t>                  &topo,
  const rd::Array2D<int>                     &label,
  const rd::Array2D<int>                     &final_label,
  const rd::Array2D<flowdir_t>               &flowdirs,
  rd::Array2D<wtd_t>                         &wtd,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  ArrayPack                                  &arp
){
 
  if(current_depression==NO_VALUE)
    return;

  auto &this_dep = deps.at(current_depression);

  //Visit child depressions. When these both overflow, then we spread water
  //across them by spreading water across their common metadepression
  MoveWaterInDepHier(this_dep.lchild, deps, topo,label,final_label,flowdirs,wtd, jump_table,arp);
  MoveWaterInDepHier(this_dep.rchild, deps, topo,label,final_label,flowdirs,wtd, jump_table,arp);

  //Catch depressions that link to the ocean through this one. These are special
  //cases because we will never spread water across the union of these
  //depressions and the current depressions: they only flow into the current
  //depression
  for(const auto c: this_dep.ocean_linked)
    MoveWaterInDepHier(c, deps, topo,label,final_label,flowdirs,wtd, jump_table,arp);

  //If the current depression is the ocean then at this point we've visited all
  //of its ocean-linked depressions (the ocean has no children). Since we do not
  //otherwise want to modify the ocean we quit here.
  if(current_depression==OCEAN)
    return;

  {
    const int lchild = this_dep.lchild;
    const int rchild = this_dep.rchild;


    //Only if both children are full should their water make its way to this
    //parent
    if(lchild!=NO_VALUE
      && deps.at(lchild).water_vol>=deps.at(lchild).wtd_vol                     
      && deps.at(rchild).water_vol>=deps.at(rchild).wtd_vol                     
    ){
    this_dep.water_vol += deps.at(lchild).water_vol + deps.at(rchild).water_vol;
    }
    assert(this_dep.water_vol >= -FP_ERROR);
    if(this_dep.water_vol < 0)
      this_dep.water_vol = 0.0;
  }


  //Each depression has an associated dep_vol. This is the TOTAL volume of the
  //meta-depression including all of its children. This property answers the
  //question, "How much water can this meta-depression hold above ground?"
  //The wtd_vol then answers how much the depression can hold in TOTAL, including below-ground. 

  //Each depression also has an associated water volume called `water_vol`. This
  //stores the TOTAL volume of the water in the depression. That is, for each
  //depression this indicates how much water is in this depression and its
  //children. However, if a depression has sufficient volume to contain its
  //water, then its water_vol will not propagate up to its parent. In this way
  //we can distinguish between depressions whose water needs to be spread versus
  //metadepressions whose children might need water spread, but which will not
  //receive any spreading themselves.

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
    this_dep.water_vol = std::min(this_dep.water_vol,this_dep.wtd_vol);
    assert(this_dep.water_vol>=-FP_ERROR);
    if(this_dep.water_vol < 0)
      this_dep.water_vol = 0.0;

    assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);
    if(this_dep.water_vol > this_dep.wtd_vol)
      this_dep.water_vol = this_dep.wtd_vol;

  } else if(this_dep.water_vol>this_dep.wtd_vol) {
    //The neighbouring depression is not the ocean and this depression is
    //overflowing (therefore, all of its children are full)
    assert(this_dep.lchild==NO_VALUE || deps.at(this_dep.lchild).water_vol==deps.at(this_dep.lchild).wtd_vol);
    assert(this_dep.rchild==NO_VALUE || deps.at(this_dep.rchild).water_vol==deps.at(this_dep.rchild).wtd_vol);

    //The marginal volume of this depression is larger than what it can hold, so
    //we determine the amount that overflows, the "extra water".
    double extra_water = this_dep.water_vol - this_dep.wtd_vol;
  
    //Now that we've figured out how much excess water there is, we fill this
    //depression to its brim. Note that we don't use addition or subtraction
    //here to ensure floating-point equality.
    this_dep.water_vol = this_dep.wtd_vol;
    assert(this_dep.water_vol>=-FP_ERROR);
    if(this_dep.water_vol < 0)
      this_dep.water_vol = 0.0;

    assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);
    if(this_dep.water_vol > this_dep.wtd_vol)
      this_dep.water_vol = this_dep.wtd_vol;


    //OverflowInto will initially send water to this depression's neighbour's
    //leaf depression via the geolink. If everything fills up, the water will
    //end up in this depression's parent. So at this point we don't have to
    //worry about the extra water here any more.

    OverflowInto(this_dep.geolink, this_dep.dep_label, this_dep.parent, topo,label,final_label,flowdirs,wtd, deps, jump_table, extra_water,arp); //TODO: use odep or geolink here?


    assert(this_dep.water_vol >= -FP_ERROR);
    if(this_dep.water_vol < 0)
      this_dep.water_vol = 0.0;

    assert(
         this_dep.water_vol==0 
      || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR
      || (this_dep.lchild==NO_VALUE && this_dep.rchild==NO_VALUE) 
      || (
              this_dep.lchild!=NO_VALUE && this_dep.rchild!=NO_VALUE
           && deps.at(this_dep.lchild).water_vol - this_dep.water_vol <= FP_ERROR
           && deps.at(this_dep.rchild).water_vol - this_dep.water_vol <= FP_ERROR
         )
    );
  }

  //All overflowing depressions should by now have overflowed all the way down
  //to the ocean. We must now spread the water in the depressions by setting
  //appropriate values for wtd.
}




///This function is similar to MoveWaterIntoPits, it is also for routing
///water physically downslope and allows infiltration as it goes, where necessary.
///When it reaches the pit cell of a depression, all the remaining water is
///moved into the DepressionHierarchy data structure for rapid flood-spill-merge
///calculations.
///Because we are only moving water in one specific depression at a time from 
///one specific inflow cell, I have used a while loop rather than a recursive call.
///
///@param extra_water   How much water, in total, is available to move downslope
///@param current_dep   The depression label we want to move water in
///@param previous_dep  The depression label the water came from
///@param topo          Topography used to generate the DepressionHierarchy
///@param wtd           Water table depth. Values of 0 indicate saturation. 
///                     Negative values indicate additional water can be added to the
///                     cell. Positive values indicate standing surface water.
///@param flowdirs      Flowdirs generated by GetDepressionHierarchy
///@param final_label   Labels from GetDepressionHierarchy indicate which top-level 
///                     (not leaf) depression each cell belongs to.
///@param label         Labels from GetDepressionHierarchy indicate which depression 
///                     each cell belongs to.
///@param deps          The DepressionHierarchy generated by GetDepressionHierarchy
///
///@return  Modifies the wtd array to reflect infiltration which has occured. 
///         Modifies the depression hierarchy `deps` to indicate the updated amount of
///         water contained in each depression. 
template<class elev_t,class wtd_t>
static void MoveWaterInOverflow(
  float                          extra_water,
  const dh_label_t               current_dep,
  const dh_label_t               previous_dep,
  const rd::Array2D<elev_t>      &topo,
  rd::Array2D<wtd_t>             &wtd,
  const rd::Array2D<flowdir_t>   &flowdirs,
  const rd::Array2D<int>         &final_label,
  const rd::Array2D<int>         &label,
  DepressionHierarchy<elev_t>    &deps,
  ArrayPack                      &arp

  ){

  auto &this_dep = deps.at(current_dep);
  auto &last_dep = deps.at(previous_dep);
  int move_to_cell;    
  int my_label = label(last_dep.out_cell);                                               //the cell that will be the next to receive the extra water
  int x,y;
  topo.iToxy(last_dep.out_cell,x,y);  

  assert(extra_water > 0);
  float infiltration_amount = extra_water*infiltration_FSM;  //Note that extra_water is a volume, so in this context infiltration_FSM is also a volume. TODO: at some locations in the code, it is a height. Standardise?
  if(extra_water <= 0.001)
    infiltration_amount = extra_water;
           
  if(infiltration_amount <= (-wtd(last_dep.out_cell))*arp.cell_area[y]){  //so the whole infiltration amount can infiltrate in this cell (but there is potentially still more water to flow downslope)
    extra_water -= infiltration_amount; 

    my_label = label(last_dep.out_cell);
      
      while(my_label != final_label(last_dep.out_cell)){ //adjust the wtd_vol of me and all my parents
        if(topo(last_dep.out_cell)<deps.at(my_label).out_elev){
          deps.at(my_label).wtd_vol -= infiltration_amount;       
          deps.at(my_label).water_vol -= infiltration_amount;
        }

        assert(deps.at(my_label).wtd_vol - deps.at(my_label).dep_vol >= -FP_ERROR);
        if(deps.at(my_label).wtd_vol < deps.at(my_label).dep_vol)
          deps.at(my_label).wtd_vol = deps.at(my_label).dep_vol;

        if(deps.at(my_label).water_vol < 0)
          deps.at(my_label).water_vol = 0; //because we might be adjusting a water_vol of a higher depression that doesn't yet have water (because we only record the water in a parent depression after the child moves its water up), so we have to switch that back to 0. 

        my_label = deps.at(my_label).parent;      

      }

      this_dep = deps.at(current_dep);

      wtd(last_dep.out_cell) += infiltration_amount/arp.cell_area[y]; //since this amount has already been subtracted from extra_water.
      //infiltration_amount was less than the available wtd space so wtd should still be negative now. 
      assert(wtd(last_dep.out_cell)<= 0);
     
  }
  else{                //Not the whole infiltration amount can infiltrate; wtd will be filled up to the level of 0 and more water will be left over. 

    extra_water += wtd(last_dep.out_cell)*arp.cell_area[y]; //+= since wtd is either negative or zero. 
   
    my_label = label(last_dep.out_cell);
      
      while(my_label != final_label(last_dep.out_cell)){ //adjust the wtd_vol of me and all my parents
        if(topo(last_dep.out_cell)<deps.at(my_label).out_elev){
          deps.at(my_label).wtd_vol += wtd(last_dep.out_cell)*arp.cell_area[y];
          deps.at(my_label).water_vol += wtd(last_dep.out_cell)*arp.cell_area[y];
        }

        assert(deps.at(my_label).wtd_vol - deps.at(my_label).dep_vol >= -FP_ERROR);
        if(deps.at(my_label).wtd_vol < deps.at(my_label).dep_vol)
          deps.at(my_label).wtd_vol = deps.at(my_label).dep_vol;

        if(deps.at(my_label).water_vol < 0)
          deps.at(my_label).water_vol = 0; //because we might be adjusting a water_vol of a higher depression that doesn't yet have water, so we have to switch that back to 0. 

        my_label = deps.at(my_label).parent;      
      }
     this_dep = deps.at(current_dep);

    wtd(last_dep.out_cell) = 0.0; //water table is now saturated in this cell. 

    //we don't need to adjust the actual infiltration_amount since that is just an indication of a portion of extra_water. We are adjusting extra_water; 
    //in the next cell we will calculate a new infiltration_amount. There is little enough space that not the whole infiltration_amount infiltrates. 
   }
    
    if(extra_water <= 0) //this would happen if above, infiltration_amount was = extra_water. 
      return;    


    int nx;
    int ny;
    if(extra_water > 0){           //Which it should be, due to if statement above, it must then move downslope but we must also check that it moves towards this_dep and doesn't flow back into last_dep
      move_to_cell = NO_VALUE;   
      for(int n=1;n<=neighbours;n++){ //Check out our neighbours
        //Use offset to get neighbour x coordinate
        nx = x+dx[n];
        //Use offset to get neighbour y coordinate
        ny = y+dy[n];      
        if(!topo.inGrid(nx,ny))  //Is cell outside grid (too far North/South)?
          continue;              //Yup: skip it.

        if((this_dep.my_subdepressions.count(label(nx,ny))!=0 || this_dep.dep_label == label(nx,ny)) && (move_to_cell == NO_VALUE || topo(nx,ny)<topo(move_to_cell)))  
          move_to_cell = topo.xyToI(nx,ny);
               
      }  //so now we know which is the correct starting cell.
           assert(move_to_cell != NO_VALUE);       
    }

  assert(extra_water > 0);

  while(extra_water > 0){
    //first let some of it evaporate for each cell it goes over. 
    extra_water = extra_water*(1-surface_evap(c));


    infiltration_amount = extra_water*infiltration_FSM;
    if(extra_water <= 0.001)
      infiltration_amount = extra_water;

    const auto ndir = flowdirs(move_to_cell);
 
    if(infiltration_amount <= -wtd(move_to_cell)*arp.cell_area[ny]){  //so the whole infiltration amount can infiltrate in this cell (but there is still more water to flow downslope)
      extra_water -= infiltration_amount;  //the whole infiltration_amount infiltrates, but there is still extra_water left over. 


      my_label = label(move_to_cell);
      
      while(my_label != final_label(move_to_cell)){ //adjust the wtd_vol of me and all my parents
        if(topo(last_dep.out_cell)<deps.at(my_label).out_elev){
          deps.at(my_label).wtd_vol -= infiltration_amount;
          deps.at(my_label).water_vol -= infiltration_amount;
        }

        assert(deps.at(my_label).wtd_vol - deps.at(my_label).dep_vol >= -FP_ERROR);
        if(deps.at(my_label).wtd_vol < deps.at(my_label).dep_vol)
          deps.at(my_label).wtd_vol = deps.at(my_label).dep_vol;

        if(deps.at(my_label).water_vol < 0)
          deps.at(my_label).water_vol = 0; //because we might be adjusting a water_vol of a higher depression that doesn't yet have water, so we have to switch that back to 0. 

        my_label = deps.at(my_label).parent;      
      }
     this_dep = deps.at(current_dep);  //I needed to reset this because seemingly the above while loops could mess it up somehow?
        
      wtd(move_to_cell) += infiltration_amount/arp.cell_area[ny];              //the cell we are in is not fully groundwater saturated. 
      assert(wtd(move_to_cell)<= FP_ERROR);
      if(wtd(move_to_cell)>=0)
        wtd(move_to_cell) = 0;
      
      if(extra_water<=0)
        break;

      assert(extra_water > 0);

      if(ndir == NO_FLOW)                   //we've reached a pit cell, so we can stop.   
        extra_water = 0;   //TODO: check that I can still just do it this way. Infiltration in this cell should still happen when I go to spread water? (check that we do infiltrate in the starting cell!)
        //TODO: confirm whether the water_vol for this depression knows about the remaining extra_water. 
        
      else{                                   //we need to keep moving downslope. 
        int x1,y1;
        topo.iToxy(move_to_cell,x1,y1);              //then we can use flowdirs to move the water all the way down to this_dep's pit cell. 
        nx = x1+dx[ndir];
        ny = y1+dy[ndir];
        move_to_cell = topo.xyToI(nx,ny);  //get the new move_to_cell
        assert(move_to_cell>=0);
      }
    }

    else{                               //there is still more infiltration_amount than will infiltrate in this cell. The wtd will become fully saturated. 
      extra_water += wtd(move_to_cell)*arp.cell_area[ny];  
      assert(extra_water > 0);
    

    my_label = label(move_to_cell);

      while(my_label != final_label(move_to_cell)){ //adjust the wtd_vol of me and all my parents
        if(topo(last_dep.out_cell)<deps.at(my_label).out_elev){
          deps.at(my_label).wtd_vol += wtd(move_to_cell)*arp.cell_area[ny];
          deps.at(my_label).water_vol += wtd(move_to_cell)*arp.cell_area[ny];
        }

        assert(deps.at(my_label).wtd_vol - deps.at(my_label).dep_vol >= -FP_ERROR);
        if(deps.at(my_label).wtd_vol < deps.at(my_label).dep_vol)
          deps.at(my_label).wtd_vol = deps.at(my_label).dep_vol;
        if(deps.at(my_label).water_vol < 0)
          deps.at(my_label).water_vol = 0; //because we might be adjusting a water_vol of a higher depression that doesn't yet have water, so we have to switch that back to 0. 

        my_label = deps.at(my_label).parent;      
        }

     this_dep = deps.at(current_dep);  //I needed to reset this because seemingly the above while loops could mess it up somehow?

      wtd(move_to_cell) = 0;              //the cell we are in is now groundwater saturated. 
      if(ndir == NO_FLOW)                   //we've reached a pit cell, so we can stop.   
        extra_water = 0;

      else{                                   //we need to keep moving downslope. 
        int x1,y1;
        topo.iToxy(move_to_cell,x1,y1);              //then we can use flowdirs to move the water all the way down to this_dep's pit cell. 
        nx = x1+dx[ndir];
        ny = y1+dy[ndir];
        move_to_cell = topo.xyToI(nx,ny);  //get the new move_to_cell
        assert(move_to_cell>=0);
      }
    }
  }
}




//When water overflows from one depression into another, this function ensures
//that chained overflows and infiltration take place; that is, that water first
//fills the water table on the path between the two depressions and then moves
//from leaf depressions to metadepressions.
//
//A depression has three places it can put the water it's received.
//The depression will try to stash water in these three places sequentially.
//  1. It can store the water in itself
//  2. It can overflow into its neighbouring depression (by following a geolink to that depression's leaf)
//  3. It can overflow into its parent
//
//Options (2) and (3) result in a recursive call. If there's enough water,
//eventually repeated calls to (3) will bring the function to the parent of the
//depression that originally called it (through its neighbour). At this point we
//stash the water in the parent and exit.
//
//Since we might end up calling the function again and again if there's a
//complex series of overflows, the `jump_table` argument holds the location of
//where the water ultimately ended up. Everything between the original node and
//this destination is then full which means that we only traverse each part of a
//filled hierarchy once.
//
//Note that since we only call this function on the leaf nodes of depressions
//the jump_table only needs to use leaves as its keys.
//
//@param root         The depression we're currently considering
//@param stop_node    When we reach this depression we dump all the excess water
//                    we're carrying into it. This depression is the parent of 
//                    the depression that first called this function. Reaching 
//                    it means that both the original metadepression and its 
//                    neighbour are both full.
///@param deps        The DepressionHierarchy generated by 
///                   GetDepressionHierarchy
///@param jump_table  A data structure that persists throughout the traversal 
///                   and is used to skip from leaf nodes to the highest known
///                   meta-depression which still has unfilled volume. Ensures 
///                   the traversal happens in O(N) time.
///@param extra_water The amount of water left to distribute. We'll try to stash
///                   it in root. If we fail, we'll pass it to root's neighbour
///                   or, if the neighbour's full, to root's parent.
//@return The depression where the water ultimately ended up. This is used to
//        update the jump table
template<class elev_t,class wtd_t>
static dh_label_t OverflowInto(
  const dh_label_t                           root,
  const dh_label_t                           previous_dep, //the previous depression. Sometimes we need to know where the water came from when we overflow it. 
  const dh_label_t                           stop_node,
  const rd::Array2D<elev_t>                  &topo,
  const rd::Array2D<int>                     &label,
  const rd::Array2D<int>                     &final_label,
  const rd::Array2D<flowdir_t>               &flowdirs,
  rd::Array2D<wtd_t>                         &wtd,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,  //Shortcut from one depression to its next known empty neighbour
  double                                     extra_water,
  ArrayPack                                  &arp 
){

  auto &this_dep = deps.at(root);
  auto &last_dep = deps.at(previous_dep);

  
  if(root==OCEAN)                        //We've reached the ocean
    return OCEAN;                        //Time to stop: there's nowhere higher in the depression hierarchy

  //FIRST PLACE TO STASH WATER: IN THIS DEPRESSION

  //We've gone around in a loop and found the original node's parent. That means
  //it's time to stop. (This may be the leaf node of another metadepression, the
  //ocean, or a standard node.)
  if(root==stop_node){                   //We've made a loop, so everything is full
    if(this_dep.parent==OCEAN){           //If our parent is the ocean
      this_dep.water_vol = this_dep.wtd_vol;
      assert(this_dep.water_vol>=-FP_ERROR);
      if(this_dep.water_vol < 0)
        this_dep.water_vol = 0;

      assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);
      return OCEAN;                      //Then the extra water just goes away
    }
    else  {                               //Otherwise
  
  //if this node is a normal parent, then I don't think it matters where in the depression the water goes, and it can just get added to this_dep.water_vol. 
  //But if the original depression is ocean_linked to this depression, then the overflow is happening from a certain location and we need to route that water.  	
  //in the case where the original depression was ocean_linked to this one, it won't be one of this depression's children. 
  	  this_dep.water_vol += extra_water; //no matter what, the water_vol needs to be updated. 
      assert(this_dep.water_vol>= -FP_ERROR);
      if(this_dep.water_vol < 0){
        this_dep.water_vol = 0.0;
      }
      if(this_dep.lchild==NO_VALUE || (this_dep.lchild != previous_dep && this_dep.rchild != previous_dep)){ //so either if this node has no children or if neither of its children was the previous depression - this implies the previous depression was ocean_linked. 
        //in this case, we should move water to the inlet of this depression and let it flow downslope. Although this is only necessary if there is groundwater space available...

        if(this_dep.wtd_vol > this_dep.dep_vol && this_dep.water_vol < this_dep.wtd_vol){ //okay, there is groundwater volume to fill, so we must move the water properly.                                                        //TODO: is the second part of that if correct?
          MoveWaterInOverflow(extra_water,this_dep.dep_label,last_dep.dep_label,topo,wtd,flowdirs,final_label,label,deps,arp);
          assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);                  
        } 

      }  
    }
    return stop_node;
  }

  if(this_dep.water_vol<this_dep.wtd_vol){                                  //Can this depression hold any water?
    const double capacity = this_dep.wtd_vol - this_dep.water_vol;          //Yes. How much can it hold?

    if(extra_water >= capacity) {                                             //It wasn't enough to hold all the water, so it will be completely filled and there is no need to worry about flow routing. 
      this_dep.water_vol = this_dep.wtd_vol;                                            //So we fill it all the way.
      assert(this_dep.water_vol >= -FP_ERROR);
      if(this_dep.water_vol < 0)
        this_dep.water_vol = 0;

      assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);

      extra_water       -= capacity;                                                    //And have that much less extra water to worry about
    }

    else{  //extra_water < capacity, so we may have to worry about routing of water and where it infiltrates. All of the extra_water will be added.  
      if((this_dep.wtd_vol > this_dep.dep_vol) && //the depression has groundwater space. We need to move water properly from the outlet. 
      ((last_dep.parent == this_dep.dep_label && last_dep.ocean_parent == true) || last_dep.parent != this_dep.dep_label)  //but only if either I'm not the last dep's parent, or I am the parent but I'm an ocean-linking parent - though these should have been handled in the above section. 
      ){                              
  
        this_dep.water_vol  = std::min(this_dep.water_vol+extra_water,this_dep.wtd_vol);      //  this_dep.water_vol += extra_water; //TODO: Is this right? Something is off about water vols and extra water but I am VERY unsure if this is right. 
        //NO - extra_water is based on what the water_vol already was!!!
        MoveWaterInOverflow(extra_water,this_dep.dep_label,last_dep.dep_label,topo,wtd,flowdirs,final_label,label,deps,arp);
        extra_water = 0;
        assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);        //TODO: Is this assert right? What assert would be appropriate here? What if thisone just needs to further overflow?
      }    

      else{ //there is either no groundwater space, or it's a parent that is not ocean-linked, so no need to worry about flow routing, just add the water to the depression. 
        this_dep.water_vol  = std::min(this_dep.water_vol+extra_water,this_dep.wtd_vol);  //Yup. But let's be careful about floating-point stuff
        assert(this_dep.water_vol>= - FP_ERROR);
        if(this_dep.water_vol < 0)
          this_dep.water_vol = 0;

        assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);
        extra_water         = 0;                                                          //No more extra water
      }
    }
  } 

  if(extra_water==0)  {                                   //If there's no more extra water
    assert(this_dep.water_vol==0 || this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);                                               
    return root;        
  }                                                                //Call it quits

  //TODO: Use jump table

  //Okay, so there's extra water and we can't fit it into this depression

  //SECOND PLACE TO STASH WATER: IN THIS DEPRESSION'S NEIGHBOUR
  //Maybe we can fit it into this depression's overflow depression!

  auto &pdep = deps.at(this_dep.parent);
  if(this_dep.odep==NO_VALUE){      //Does the depression even have such a neighbour? 

    if(this_dep.parent!=OCEAN && pdep.water_vol==0) //At this point we're full and heading to our parent, so it needs to know that it contains our water
      pdep.water_vol += this_dep.water_vol;
    this_dep.water_vol = this_dep.wtd_vol;
    return jump_table[root] = OverflowInto(this_dep.parent, this_dep.dep_label, stop_node, topo,label,final_label,flowdirs,wtd, deps, jump_table, extra_water,arp);  //Nope. Pass the water to the parent
  }

  //Can overflow depression hold more water?
  auto &odep = deps.at(this_dep.odep);
  if(odep.water_vol<odep.wtd_vol){  //Yes. Move the water geographically into that depression's leaf.

    if(this_dep.parent!=OCEAN && pdep.water_vol==0 && odep.water_vol+extra_water>odep.wtd_vol) //It might take a while, but our neighbour will overflow, so our parent needs to know about our water volumes
      pdep.water_vol += this_dep.water_vol + odep.wtd_vol;           //Neighbour's water_vol will equal its dep_vol
  //TODO: is this right? Aren't we adding the odep.wtd_vol here AND then ading it again when we actually overflow into the parent?
    deps.at(this_dep.geolink).water_vol += extra_water; //TODO: and should I also subtract the extra water from this_dep?
    this_dep.water_vol = this_dep.wtd_vol;

    //TODO: No I shouldn't, since it gets set equal to wtd vol before I call overflowinto. However, do I need to sometimes subtract it? Subtract it somewhere else?
   // this_dep.water_vol -= extra_water;
   // this_dep.water_vol = std::max(this_dep.water_vol,0.0);
    return jump_table[root] = OverflowInto(this_dep.geolink, this_dep.dep_label, stop_node, topo,label,final_label,flowdirs,wtd, deps, jump_table, extra_water,arp);
  }  //TODO: I am concerned that using the geolink here may actually no longer be the best choice now that I've implemented downslope flow of water when overflowing. 
     //Imagine the case where dep A is geolinked to dep B but dep B is a part of metadep C and A and B don't actually touch. During overflow, we 
     //won't find any cells of B to flow into. So I think linking to C, i.e. the odep, will actually work best?
     //Geolinks were needed when we were just tossing the water into the bottom of the depression but seem wrong with the new functionality. 

  //Okay, so the extra water didn't fit into this depression or its overflow
  //depression. That means we pass it to this depression's parent.

  //If we've got here we have a neighbour, but we couldn't stash water in the
  //neighbour because it was full. So we need to see if our parent knows about
  //us.
  if(this_dep.parent!=OCEAN && pdep.water_vol==0)
    pdep.water_vol += this_dep.water_vol + odep.water_vol;
    
  this_dep.water_vol = this_dep.wtd_vol;

  //THIRD PLACE TO STASH WATER: IN THIS DEPRESSION'S PARENT
  return jump_table[root] = OverflowInto(this_dep.parent, this_dep.dep_label, stop_node,topo,label,final_label,flowdirs,wtd, deps, jump_table, extra_water,arp);
}



//When we're determining which depressions to spread standing surface water
//across, we need to know a leaf depression to start filling from and the
//metadepression that contains all of the water. We also need to know all the
//depressions across which we're allowed to spread water. This data structure
//keeps track of this information.
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
};



///This function recursively traverses the depression hierarchy until it reaches
///a leaf depression. It then starts climbing back up. If a leaf depression is
///full, it notes the leaf depression's id. This a potential place to start
///trying to flood a metadepression.
///
///At higher levels if both of a metadepression's children are full, then the
///metadepression makes a note of the ids of all the depressions contained
///within their subtrees. It also chooses arbitrarily one of the leaf
///depressions to start flooding from (since both children are full, the choice
///of a starting point for flooding doesn't matter).
///
///Eventually, a metadepression is reached which has more volume than water. At
///this point a helper function `FillDepressions()` is called. This is passed
///the ids of the arbitrarily chosen leaf depression, the partially-filled
///metadepression, and all the depressions in between. It then modifies the
///water table so that standing water rises to its natural level within the
///partially-filled metadepression.
///
///@param current_depression  The depression we're currently considering
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression 
///                each cell belongs to.
///@param wtd      Water table depth. Values of 0 indicate saturation. 
///                Negative values indicate additional water can be added to the
///                cell. Positive values indicate standing surface water.
///@return Information about the subtree: its leaf node, depressions it 
///        contains, and its root node.
template<class elev_t>
static SubtreeDepressionInfo FindDepressionsToFill(
  const int                         current_depression,    //Depression we are currently in
  const DepressionHierarchy<elev_t> &deps,                  //Depression hierarchy
  const rd::Array2D<float>          &topo,                  //Topographic data (used for determinining volumes as we're spreading stuff)
  const rd::Array2D<dh_label_t>     &label,                 //Array indicating which leaf depressions each cell belongs to
  rd::Array2D<float>                &wtd,                    //Water table depth
  ArrayPack                         &arp
){
  //Stop when we reach one level below the leaves
  if(current_depression==NO_VALUE)
    return SubtreeDepressionInfo();

  const auto& this_dep = deps.at(current_depression);

  //We start by visiting all of the ocean-linked depressions. They don't need to
  //pass us anything because their water has already been transferred to this
  //metadepression tree by MoveWaterInDepHier(). Similarly, it doesn't matter what their leaf
  //labels are since we will never spread water into them.
  for(const auto c: this_dep.ocean_linked)
    FindDepressionsToFill(c, deps, topo, label, wtd,arp);

  //At this point we've visited all of the ocean-linked depressions. Since all
  //depressions link to the ocean and the ocean has no children, this means we
  //have visited all the depressions and spread their water. Since we don't wish
  //to modify the ocean, we are done.
  if(current_depression==OCEAN)
    return SubtreeDepressionInfo();

  //We visit both of the children. We need to keep track of info from these
  //because we may spread water across them.
  SubtreeDepressionInfo left_info  = FindDepressionsToFill(this_dep.lchild, deps, topo, label, wtd,arp);
  SubtreeDepressionInfo right_info = FindDepressionsToFill(this_dep.rchild, deps, topo, label, wtd,arp );   

  SubtreeDepressionInfo combined;
  combined.my_labels.emplace(current_depression);
  combined.my_labels.merge(left_info.my_labels);
  combined.my_labels.merge(right_info.my_labels);

  combined.leaf_label = left_info.leaf_label;  //Choose left because right is not guaranteed to exist
  if(combined.leaf_label==NO_VALUE)            //If there's no label, then there was no child
    combined.leaf_label = current_depression;  //Therefore, this is a leaf depression

  combined.top_label = current_depression;

  //The water volume should never be greater than the depression volume because
  //otherwise we would have overflowed the water into the neighbouring
  //depression and moved the excess to the parent.
  assert(this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);
 // if(this_dep.water_vol > this_dep.wtd_vol)
   // this_dep.water_vol = this_dep.wtd_vol; 

  //Since depressions store their marginal water volumes, if a parent depression
  //has 0 marginal water volume, then both of its children have sufficient
  //depression volume to store all of their water. TODO: is this always true? What about a metadepression that is the EXACT size of its two children, so it never has a marginal water volume?
  // However, if our parent is an
  //ocean-link then we are guaranteed to be able to fill now because excess
  //water will have been transferred into the parent and we don't want to pool
  //the parent's water with our own (it might be at the bottom of a cliff).

//  if(this_dep.water_vol==this_dep.dep_vol && deps.at(this_dep.parent).water_vol==0){
  //  std::cout<<"here is an exact equal depression"<<this_dep.dep_label<<std::endl;
  //  for(int y=0;y<label.height();y++)
    //  for(int x=0;x<label.width();x++){
      //  if(label(x,y)==this_dep.dep_label){   //TODO: Will this catch everything that is part of this depression?
        //  if(topo(x,y)<this_dep.out_elev){
          //  assert(wtd(x,y)==0);
       //     wtd(x,y) = this_dep.out_elev - topo(x,y);
      //    }
      //    }
      //  }
   // return SubtreeDepressionInfo();

   // }


  if(this_dep.water_vol<this_dep.wtd_vol || this_dep.ocean_parent || (this_dep.water_vol==this_dep.dep_vol && deps.at(this_dep.parent).water_vol==0)){  //TODO: Should there also be an OR if this dep's parent is the ocean?
    assert(this_dep.water_vol - this_dep.wtd_vol <= FP_ERROR);


 //   if(this_dep.water_vol > this_dep.wtd_vol)
   //   this_dep.water_vol = this_dep.wtd_vol;
  //  assert(deps.at(this_dep.parent).dep_label == OCEAN || deps.at(this_dep.parent).water_vol == 0);
     //I no longer understand the above assert. Why on earth would the parent of an ocean-linked depression not be allowed to have a non-zero water vol?


    //If both of a depression's children have already spread their water, we do not
    //want to attempt to do so again in an empty parent depression. 
    //We check to see if both children have finished spreading water. 
  
     FillDepressions(combined, this_dep.water_vol, deps, topo, label, wtd,arp);

    //At this point there should be no more water all the way up the tree until
    //we pass through an ocean link, so we pass this up as a kind of null value.
    return SubtreeDepressionInfo();
  } else {
    return combined;
  }
}



///This function adjusts the water table to reflect standing surface water that
///has pooled at the bottom of depressions.
///
///At this point we know the id of an arbitrarily-chosen leaf depression and the
///id of the metadepression across which we will spread the water.
///
///Starting at the pit cell of the leaf depression we use a priority queue to
///climb from the lowest cells to the highest cell. For each cell we determine
///how much volume is contained within the depression if it were flooded to the
///level of this cell. If the volume is sufficent, we do the flooding.
///
///Therefore, we need an efficient way to determine the volume of a depression.
///To do so, we note that if the elevation of the current cell is $o$ and a
///depression contains cells of elevations $\{a,b,c,d\}$, then the volume of the
///depression is $(o-a)+(o-b)+(o-c)+(o-d)=4o-a-b-c-d=4o-sum(elevations)$. Thus,
///if we keep track of the number of cells in a depression and their total
///elevation, it is possible to calculate the volume of a depression at any time
///based on a hypothetical outlet level. We call this the Water Level Equation.
///
///@param stdi     Leaf node, metadepression, depressions between. Determines 
///                the extent of the flooding.
///@param water_vol How much water needs to be spread throughout the depression
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression 
///                each cell belongs to.
///@param wtd      Water table depth. Values of 0 indicate saturation. 
///                Negative values indicate additional water can be added to the
///                cell. Positive values indicate standing surface water.
///@return         N/A
template<class elev_t, class wtd_t>
static void FillDepressions(
  //Identifies a meta-depression through which water should be spread, leaf node
  //from which the water should be spread, and valid depressions across which
  //water can spread.
  SubtreeDepressionInfo             &stdi,  
  double                            water_vol, //Amount of water to spread around this depression
  const DepressionHierarchy<elev_t> &deps,      //Depression hierarchy
  const rd::Array2D<float>          &topo,      //Topographic data for calculating marginal volumes as we attempt to spread water
  const rd::Array2D<dh_label_t>     &label,     //2D array in which each cell is labeled with the leaf depression it belongs to
  rd::Array2D<wtd_t>                &wtd,        //Water table depth: we transfer water into this
  ArrayPack                         &arp
){
  //Nothing to do if we have no water
  if(water_vol==0)
    return;


  //Hashset stores the ids of the cells we've visited. We don't want to use a 2D
  //array because the size of the DEM as a whole could be massive and that's a
  //lot of memory to allocate/deallocate each time this function is called. We
  //arbitrarily reserve enough space in the hashset for a few thousand items.
  //This should be large than most depressions while still being small by the
  //computer's standards.

  std::unordered_set<int> visited(2048);

  //Priority queue that sorts cells by lowest elevation first. If two cells are
  //of equal elevation the one added most recently is popped first. The ordering
  //of the cells processed by this priority queue need not match the ordering of
  //the cells processed by the depression hierarchy.
  rd::GridCellZk_high_pq<elev_t> flood_q;  
  rd::GridCellZk_high_pq<elev_t> neighbour_q;    
  
  { //Scope to limit pit_cell
    //Cell from which we can begin flooding the meta-depression. Which one we
    //choose is arbitrary, since we will fill all of the leaf depressions and
    //intermediate metadepressions until and including when we reach the
    //depression identified by stdi.top_label
    const auto pit_cell    = deps.at(stdi.leaf_label).pit_cell;
    assert(pit_cell>=0);

    //We start flooding at the pit cell of the depression and work our way
    //upwards
    flood_q.emplace(
      pit_cell % topo.width(),
      pit_cell / topo.width(),
      topo(pit_cell)
    );

    visited.emplace(pit_cell);
  }

  //Cells whose wtd will be affected as we spread water around
  std::vector<int> cells_affected;

  //Stores the sum of the elevations of all of the cells in cells_affected. Used
  //for calculating the volume we've seen so far. (See explanation above or in
  //dephier.hpp TODO)
  double current_elevation = 0;
  double previous_elevation = 0;

//Stores the sum of the elevations of all of the cells in cells_affected. Used
  //for calculating the volume we've seen so far. (See explanation above or in
  //dephier.hpp TODO)
  double total_elevation = 0;

  double current_volume = 0;
  double current_area = 0;


  while(!flood_q.empty()){
    const auto c = flood_q.top();
    flood_q.pop();
    current_elevation = static_cast<double>(topo(c.x,c.y));

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
    //current_volume = cells_affected.size()*static_cast<double>(topo(c.x,c.y)) - total_elevation; 

    current_volume += (current_elevation-previous_elevation)*current_area;

    assert(water_vol >= - FP_ERROR);
    if(water_vol < 0)
      water_vol = 0; 

    //All the cells within this depression should have water table depths less
    //than or equal to zero because we have moved all of their water down slope
    //into the pit cell. Since we may already have filled other depressions
    //their cells are allowed to have wtd>0. Thus, we raise a warning if we are
    //looking at a cell in this unfilled depression with wtd>0.
    if(stdi.my_labels.count(label(c.x,c.y))==1 && wtd(c.x,c.y)>0)
      throw std::runtime_error("A cell was discovered in an unfilled depression with wtd>0!");

    //There are two possibilities:
    //1. The virtual water level exceeds slightly the height of the cell. The cell's water table then fills up as much as it can.
    //   The water surface is then level with the height of the cell.
    //
    //2. The water surface is below the height of the cell because there is sufficient topographic volume to hold all the water.
    //   In this case, the cell's water table is left unaffected.

    if(water_vol<=current_volume-wtd(c.x,c.y)*arp.cell_area[c.y]){
      //The current scope of the depression plus the water storage capacity of
      //this cell is sufficient to store all of the water. We'll stop adding
      //cells and fill things now.

      //We will fill the depression so that the surface of the water is at this
      //elevation.
      double water_level;

      if(current_volume<water_vol){ //TODO: Check stdi.my_labels.count(label(c.x,c.y))==0 ?
        //The volume of water exceeds what we can hold above ground, so we will
        //stash as much as we can in this cell's water table. This is okay
        //because the above ground volume plus this cell's water table IS enough
        //volume (per the if-clause above).


        //Fill in as much of this cell's water table as we can
        const double fill_amount = water_vol - current_volume;                               //TODO is it needed to update wtd_vol here? I think the code can work without doing so, but it would be better for error checking if we did. If we do, best way to do so?
        assert(fill_amount >= -FP_ERROR);
   //     if(fill_amount < 0)
     //     fill_amount = 0; 

        wtd(c.x,c.y)   += fill_amount/arp.cell_area[c.y];
        water_vol -= fill_amount;   //Doesn't matter because we don't use water_vol anymore
        water_level     = topo(c.x,c.y);

      } else if (current_volume==water_vol)   //The volume of water is exactly equal to the above ground volume so we set the water level equal to this cell's elevation
          water_level = topo(c.x,c.y);
      else {  //The water volume is less than this cell's elevation, so we calculate what the water level should be.
        //We have that Volume = (Water Level)*(Cell Count)-(Total Elevation)
        //rearranging this gives us:
   //     water_level = (water_vol+total_elevation)/cells_affected.size();//TODO
   //       water_level = (water_vol/current_area)+total_elevation/cells_affected.size();
          water_level = topo(c.x,c.y);  //TODO: THIS IS INCORRECT! The real value should be something less than topo(c.x,c.y) but more than the topo of the previous cell. Not sure how to calculate it now that water_vol is an actual volume...

      }
      //Water level must be higher than (or equal to) the previous cell we looked at, but lower than (or equal to) the current cell
     // std::cout<<"topo "<<current_elevation<<" "<<previous_elevation<<std::endl;
     // std::cout<<" water level "<<water_level<<std::endl;
     // std::cout<<"current area "<<current_area<<" water vol "<<water_vol<<" total elevation "<<total_elevation<<" cells affected size "<<cells_affected.size()<<std::endl;
      assert(cells_affected.size()==0 || topo(cells_affected.back()) - water_level <= FP_ERROR); 
      assert(topo(c.x,c.y)-water_level >= -FP_ERROR);

      for(const auto c: cells_affected){
        assert(wtd(c) >= -FP_ERROR);               //This should be true since we have been filling wtds as we go.
        if(wtd(c) < 0)
          wtd(c) = 0;

          
        assert(water_level - topo(c) >= -FP_ERROR);
        if(water_level < topo(c))
          water_level = topo(c);
        
        wtd(c) = water_level - topo(c);  //only change the wtd if it is an increase, here. We can't take water away from cells that already have it (ie reduce groundwater in saddle cells within a metadepression.)
        if(-FP_ERROR<=wtd(c) && wtd(c)<0)
          wtd(c) = 0;
        assert(wtd(c)>= -FP_ERROR);
        if(wtd(c) < 0)
          wtd(c) = 0;
      }
      //We've spread the water, so we're done        
      return;
      
    }  else {
      //We haven't found enough volume for the water yet.

      //During the adding of neighbours neighbours might get added that are
      //lower than we are and belong to a different depression (notably, this
      //happens at the edge of a flat abuting an ocean). These cells will then
      //be popped and could be processed inappropriately. To prevent this, we
      //skip them here.
      if(stdi.my_labels.count(label(c.x,c.y))==0){  //CHECK. This was preventing cells that flowed to the ocean from allowing my depression volume to update. Is this way ok? Is this even needed?
        continue;
      }
    
      //Okay, we're allow to add this cell's neighbours since this cell is part
      //of the metadepression.

      //Add this cell to those affected so that its volume is available for
      //filling.
      cells_affected.emplace_back(topo.xyToI(c.x,c.y));
      //Fill in cells' water tables as we go

      assert(wtd(c.x,c.y) <= FP_ERROR);
      if(wtd(c.x,c.y) > 0)
        wtd(c.x,c.y) = 0;
      water_vol += wtd(c.x,c.y)*arp.cell_area[c.y];  //We use += because wtd is less than or equal to zero
      wtd(c.x,c.y)    = 0;             //Now we are sure that wtd is 0, since we've just filled it

           //Add the current cell's information to the running total
      total_elevation += topo(c.x,c.y);
      current_area += arp.cell_area[c.y];  //adding to the area after the volume because when there is only 1 cell, the answer for volume should be 0, etc. Don't want to include the area of the target cell. 


      
 
      for(int n=1;n<=neighbours;n++){
        const int nx = c.x + dx[n]; //TODO ModFloor(x+dx[n],topo.width()); //Get neighbour's x-coordinate using an offset and wrapping
        const int ny = c.y + dy[n];                     //Get neighbour's y-coordinate using an offset
        if(!topo.inGrid(nx,ny))                         //Is this cell in the grid?
          continue;                                     //Nope: out of bounds.
        const int ni = topo.xyToI(nx,ny);               //Get neighbour's flat index
     
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
  
        if(visited.count(ni)==0){ //&& ((label(nx,ny)!=OCEAN) || arp.land_mask(nx,ny)>0.0f)){
          if(stdi.my_labels.count(label(nx,ny))==0)  //CHECK. This was preventing cells that flowed to the ocean from allowing my depression volume to update. Is this way ok? Is this even needed?
            neighbour_q.emplace(nx,ny,topo(nx,ny));
          else
            flood_q.emplace(nx,ny,topo(nx,ny));
          visited.emplace(ni);
        }
      }
    }


    if(flood_q.empty() && !neighbour_q.empty()){
      const auto c = neighbour_q.top();
      neighbour_q.pop();    
      flood_q.emplace(c.x,c.y,topo(c.x,c.y));
    }

    previous_elevation = static_cast<double>(topo(c.x,c.y));
  }

  //Since we're in this function we are supposed to be guaranteed to be able to
  //fill our depression, since we only enter this function once that is true.
  //Therefore, if we've reached this point, something has gone horribly wrong
  //somewhere. :-(

  std::cerr<<"E PQ loop exited without filling a depression!"<<std::endl;
  std::cerr<<"water vol "<<water_vol<<" current vol "<<current_volume<<std::endl;

  assert(!flood_q.empty());
  throw std::runtime_error("PQ loop exited without filling a depression!");
}

}

#endif
