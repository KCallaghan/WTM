#pragma once

#include "dephier.hpp"
#include "parameters.hpp"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/math.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/timer.hpp>

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
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace richdem::dephier {

constexpr auto dx = d8x;
constexpr auto dy = d8y;
constexpr auto dinverse = d8_inverse;
constexpr int        neighbours = 8;

constexpr double FP_ERROR = 1e-4;

///////////////////////////////////
//Function prototypes
///////////////////////////////////

template<class elev_t>
void ResetDH(DepressionHierarchy<elev_t> &deps);

template<class elev_t>
void FillSpillMerge(
  Parameters                    &params,
  DepressionHierarchy<elev_t>   &deps,
  ArrayPack                     &arp
);

template<class elev_t>
static void MoveWaterIntoPits(
  Parameters                    &params,
  DepressionHierarchy<elev_t>   &deps,
  ArrayPack                     &arp
);

template<class elev_t, class wtd_t>
static void CalculateWtdVol(
  rd::Array2D<wtd_t>            &wtd,
  const std::vector<double>     &cell_area,
  const rd::Array2D<float>      &porosity,
  const rd::Array2D<dh_label_t> &final_label,
  DepressionHierarchy<elev_t>   &deps
);

static double CalculateInfiltration(
  int                           cell,
  double                        distance,
  double                        h_0,
  ArrayPack                     &arp
);

template<class elev_t>
static void MoveWaterInDepHier(
  int                                        current_depression,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  Parameters                                 &params,
  ArrayPack                                  &arp
);

template<class elev_t>
static void MoveWaterInOverflow(
  double                         extra_water,
  const dh_label_t               current_dep,
  const dh_label_t               previous_dep,
  DepressionHierarchy<elev_t>    &deps,
  Parameters                     &params,
  ArrayPack                      &arp
  );

template<class elev_t>
static dh_label_t OverflowInto(
  const dh_label_t                           root,
  const dh_label_t                           previous_dep,
  const dh_label_t                           stop_node,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  double                                     extra_water,
  Parameters                                 &params,
  ArrayPack                                  &arp
);

class SubtreeDepressionInfo;

template<class elev_t>
static SubtreeDepressionInfo FindDepressionsToFill(
  const int                         current_depression,
  DepressionHierarchy<elev_t>       &deps,
  ArrayPack                         &arp
);

template<class elev_t>
static void FillAFullDepression(
    SubtreeDepressionInfo             &stdi,
    DepressionHierarchy<elev_t>       &deps,
    ArrayPack                         &arp
);

template<class elev_t>
static void FillDepressions(
  const int                            pit_cell,
  const int                            out_cell,
  const double                         out_elev,
  const std::unordered_set<dh_label_t> &dep_labels,
  double                               water_vol,
  DepressionHierarchy<elev_t>          &deps,
  ArrayPack                            &arp
);

template<class elev_t, class wtd_t>
void BackfillDepression(
  double water_level,
  const Array2D<elev_t>   &topo,
  Array2D<wtd_t>          &wtd,
  std::vector<flat_c_idx> &cells_affected
);

template<class elev_t, class wtd_t>
double DetermineWaterLevel(
  const Array2D<elev_t>     &topo,
  const Array2D<float>      &porosity,
  const std::vector<double> &cell_area,
  Array2D<wtd_t>            &wtd,
  const double water_vol,
  const int    cx,
  const int    cy,
  const double current_volume,
  const double current_area,
  const double area_times_elevation_total
);


///////////////////////////////////
//Implementations
///////////////////////////////////



///This function routes surface water into pit cells and then distributes
///it so that it fills the bottoms of depressions, taking account of overflows.
///
///@param params   Various paramaters that are used here and in other portions
///                of the code. In FSM, we use:
///                - params.infiltration_on: True or False for if the user
///                  wants infiltration to occur as water flows downslope.
///                - params.cellsize_n_s_metres: cellsize in the north-south
///                  direction, used to calculate a distance travelled by
///                  flowing water when infiltration is turned on.
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param arp      ArrayPack consisting of all of the various arrays used by this
///                and the other portions of the code. In FSM, we use:
///                - arp.cell_area: Area of a cell, important for calculating
///                  the total volume of water held in that cell.
///                - arp.cellsize_e_w_metres: cellsize in the east-west direction,
///                  used to calculated a distance travelled by flowing water
///                  when infiltration is turned on.
///                - arp.final_label: the uppermost label (highest metadepression)
///                  to which a cell belongs.
///                - arp.flowdirs: the flow directions calculated in DH.
///                - arp.infiltration_array: records the amount of infiltration
///                  that has occurred in total during downslope flow, if
///                  infiltration is turned on, for informational purposes.
///                - arp.label: the label of the leaf depression associated
///                  with each cell.
///                - arp.porosity: the porosity associated with each cell.
///                  We assume that porosity does not change with depth.
///                - arp.runoff: used to record the amount of surface water in
///                  each cell, when infiltration is turned on. This allows us
///                  to have only a portion of water infiltrate, and still have
///                  surface water present.
///                - arp.slope: slope values for each cell. Used to calculate
///                  velocity of water flow when infiltration is turned on.
///                - arp.topo: the topography, i.e. elevation above sea level
///                  of each cell.
///                - arp.vert_ksat: the vertical hydraulic conductivity of
///                  each cell, used in calculations of infiltration when
///                  this is turned on.
///                - arp.wtd: the depth to water table in each cell. This is the
///                  array which is updated by FSM once we have calculated how
///                  and where water moves, and which depressions contain
///                  surface water.
///
///Note that from GetDepressionHierarchy we already know the number of cells
///within each depression as well as its volume.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water contained in each depression. Modifies `wtd` to indicate how
///        saturated a cell is or how much standing surface water it has.
template<class elev_t>
void FillSpillMerge(
  Parameters                    &params,
  DepressionHierarchy<elev_t>   &deps,
  ArrayPack                     &arp
){
  Timer timer_overall;
  timer_overall.start();

  RDLOG_ALG_NAME<<"FillSpillMerge";

  //If we're using the DH more than once we need to reset its water variables
  //between uses.
  ResetDH(deps);

  //We move standing water downhill to the pit cells of each depression
  MoveWaterIntoPits(params, deps, arp);

  //Calculate the wtd_vol of depressions, in order to be able to know which
  //need to overflow and which can accommodate more water:
  //the reason for only calculating it here is that it's better to do this
  //after MoveWaterIntoPits, when a lot of the infiltration
  //that would occur has already happened (if the user has infiltration turned on),
  //so we didn't have to constantly update wtd_vol during MoveWaterIntoPits.
  //wtd_vol is the total water that a depression can accommodate, including the
  //above ground depression volume plus any below-ground groundwater space.
  CalculateWtdVol(arp.wtd, arp.cell_area, arp.porosity, arp.final_label, deps);

  {
    //Scope to limit `timer_overflow` and `jump_table`. Also ensures
    //`jump_table` frees its memory
    Timer timer_overflow;
    timer_overflow.start();
    std::unordered_map<dh_label_t, dh_label_t> jump_table;
    //Now that the water is in the pit cells, we move it around so that
    //depressions which contain too much water overflow into depressions that
    //have less water. If enough overflow happens, then the water is ultimately
    //routed to the ocean.
    MoveWaterInDepHier(OCEAN, deps,jump_table,params,arp);
    RDLOG_TIME_USE<<"FlowInDepressionHierarchy: Overflow time = "<<timer_overflow.stop();
  }

  //Sanity checks
  for(int d=1;d<(int)deps.size();d++){
    const auto &dep = deps.at(d);

    assert(dep.water_vol==0 || fp_le(dep.water_vol,dep.wtd_vol) );
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.lchild!=NO_VALUE && fp_le(deps.at(dep.lchild).water_vol,dep.water_vol) ));
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.rchild!=NO_VALUE && fp_le(deps.at(dep.rchild).water_vol,dep.water_vol) ));

    if(dep.lchild != NO_VALUE && deps.at(dep.lchild).water_vol > dep.water_vol && dep.water_vol>0)//fp_ge(dep.water_vol,0))
      deps.at(dep.lchild).water_vol = dep.water_vol;
    if(dep.rchild != NO_VALUE && deps.at(dep.rchild).water_vol > dep.water_vol  && dep.water_vol>0)//fp_ge(dep.water_vol,0))
      deps.at(dep.rchild).water_vol = dep.water_vol;
  }

  RDLOG_PROGRESS<<"p Finding filled...";
  rd::Timer timer_filled;
  timer_filled.start();
  //We start at the ocean, crawl to the bottom of the depression hierarchy and
  //determine which depressions or metadepressions contain standing water. We
  //then modify `wtd` in order to distribute this water across the cells of the
  //depression which will lie below its surface.
  std::cout<<"find depressions"<<std::endl;
  FindDepressionsToFill(OCEAN,deps,arp);
  RDLOG_TIME_USE<<"t FlowInDepressionHierarchy: Fill time = "<<timer_filled.stop()<<" s";
  RDLOG_TIME_USE<<"t FlowInDepressionHierarchy = "<<timer_overall.stop()<<" s";
}



///This function works much like a standard flow accumulation algorithm except
///that as the water moves downhill it contributes to saturating the water table
///and, when it reaches the pit cell of a depression, all the remaining water is
///moved into the DepressionHierarchy data structure for rapid flood-spill-merge
///calculations.
///If infiltration is turned on, we calculate a proportion of water to infiltrate
///based on its velocity moving across the cell, and on how far it must travel.
///If infiltration is turned off, we skip the cell-by-cell movement of water,
///and move the water directly into the DepressionHierarchy.
///
///@param params   Various paramaters that are used here and in other portions
///                of the code. Here we use:
///                - params.infiltration_on, to determine whether infiltration
///                  and cell-by-cell water movement should occur, and
///                - params.cellsize_n_s_metres, to determine the distance
///                  travelled by water, if infiltration is on.
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param arp      ArrayPack consisting of all of the various arrays used by this
///                and the other portions of the code. Here we use:
///                - arp.cell_area: Area of a cell, important for calculating
///                  the total volume of water held in that cell.
///                - arp.cellsize_e_w_metres: cellsize in the east-west direction,
///                  used to calculated a distance travelled by flowing water
///                  when infiltration is turned on.
///                - arp.flowdirs: the flow directions calculated in DH.
///                - arp.infiltration_array: records the amount of infiltration
///                  that has occurred in total during downslope flow, if
///                  infiltration is turned on, for informational purposes.
///                - arp.label: the label of the leaf depression associated
///                  with each cell.
///                - arp.runoff: used to record the amount of surface water in
///                  each cell, when infiltration is turned on. This allows us
///                  to have only a portion of water infiltrate, and still have
///                  surface water present.
///                - arp.topo: the topography, i.e. elevation above sea level
///                  of each cell.
///                - arp.wtd: the depth to water table in each cell. If the water
///                  table is below the land surface (i.e. groundwater), then
///                  wtd will be negative. Positive values indicate the
///                  presence of surface water.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water contained in each LEAF depression. Modifies `wtd` to indicate
///        how saturated a cell is. All values in wtd will be <=0 following this
///        operation.
template<class elev_t>
static void MoveWaterIntoPits(
  Parameters                   &params,
  DepressionHierarchy<elev_t>  &deps,
  ArrayPack                    &arp
){
  rd::Timer timer;
  rd::ProgressBar progress;
  timer.start();

  RDLOG_PROGRESS<<"Moving surface water downstream...";

  //if infiltration is turned off, then we can skip most of the work in this function!
  //there is no need to trace the cell-to-cell path of water downslope if
  //infiltration will not occur. Since we already know the labels of each cell, we
  //know which depression the water belongs in.
  //So, if infiltration is turned off, we can just move this water directly
  //into the depression hierarchy.
  if(params.infiltration_on == false){
    //move the water to the appropriate depression's water_vol.
    for(int y=0;y<arp.topo.height();y++)
    for(int x=0;x<arp.topo.width(); x++){
      if(arp.wtd(x,y)>0){             //There is surface water in the cell.
        if(arp.label(x,y) == OCEAN){  //If downstream neighbour is the ocean, we drop our water into it and the ocean is unaffected.
          deps.at(OCEAN).water_vol += arp.wtd(x,y)*arp.cell_area[y]; //recording the volume that has ended up in the ocean.
          params.total_loss_to_ocean += arp.wtd(x,y)*arp.cell_area[y];
          arp.wtd(x,y) = 0.0;         //all surface water from the cell has flowed into the ocean.
        }
        else{
          deps[arp.label(x,y)].water_vol += arp.wtd(x,y)*arp.cell_area[y];
          //wtd is a depth of water and water_vol is a volume, so multiply by the area of the cell to convert.
          assert(fp_ge(deps[arp.label(x,y)].water_vol,0));
          if(deps[arp.label(x,y)].water_vol < 0)
            deps[arp.label(x,y)].water_vol = 0.0;
          arp.wtd(x,y) = 0; //Clean up as we go: all surface water from the cell has now been recorded in the depression hierarchy.
        }
      }
    }
  }

  //If infiltration is turned on, then we need to move all of the water
  //downstream into pit cells. To do so, we use the steepest-descent
  //flow directions provided by the depression hierarchy code.

  //First, move any surface water into a separate array.
  //This is to allow us to have infiltration into wtd, that does not
  //necessarily fill all the way to the surface,
  //and still have water continue moving downslope.
  else{   //infiltration is turned on, so we have to move water cell by cell and calculate how much will infiltrate as we go.
    #pragma omp parallel for default(none) shared(arp)
    for(unsigned int i=0;i<arp.topo.size();i++){
      if(arp.wtd(i)>0){
        arp.runoff(i) = arp.wtd(i);
        arp.wtd(i) = 0;                //And we reset any cells that contained surface water to 0: this is now stored in arp.runoff.
        arp.infiltration_array(i) = 0; //The infiltration array is informational, to record the total infiltration that has occured
                                       //in one round of FSM, so here we reset this to 0.
      }
    }

    //Calculate how many upstream cells flow into each cell
    rd::Array2D<char>  dependencies(arp.topo.width(),arp.topo.height(),0);
    #pragma omp parallel for default(none) shared(arp, dependencies, dinverse, dx, dy) collapse(2)
    for(int y=0;y<arp.topo.height();y++)
    for(int x=0;x<arp.topo.width(); x++)
    for(int n=1;n<=neighbours;n++){         //Loop through neighbours
      const int nx = x+dx[n];               //Identify coordinates of neighbour
      const int ny = y+dy[n];
      if(!arp.topo.inGrid(nx,ny))           //TODO: add east-west wraparound
        continue;
      if(arp.flowdirs(nx,ny)==dinverse[n])  //Does my neighbour flow into me?
        dependencies(x,y)++;                //Increment my dependencies
    }

    //Find the peaks. These are the cells into which no other cells pass flow
    //(i.e. 0 dependencies). We know the flow accumulation of the peaks
    //without having to perform any recursive calculations;
    //they just pass flow downstream. From the peaks, we
    //can begin a breadth-first traversal in the downstream direction by adding
    //each cell to the frontier/queue as its dependency count drops to zero.
    std::queue<flat_c_idx> q;
    for(unsigned int i=0;i<arp.topo.size();i++){
      if(dependencies(i)==0)                //Is it a peak?
        q.emplace(i);
    }  //Yes.


    double distance = 0;  //this is used below when calculating how much infiltration will occur.
    //Starting with the peaks, pass flow downstream
    progress.start(arp.topo.size());
    while(!q.empty()){
      ++progress;

      const auto c = q.front();          //Copy focal cell from queue
      q.pop();                           //Clear focal cell from queue

      //Coordinates of downstream neighbour, if any
      const auto ndir = arp.flowdirs(c); //Neighbour direction
      int n = NO_FLOW;                   //Neighbour address
      int x, y, nx, ny;
      arp.topo.iToxy(c,x,y);
      if(ndir!=NO_FLOW){
        nx = x+dx[ndir];
        ny = y+dy[ndir];
        n  = arp.topo.xyToI(nx,ny);
        assert(n>=0);
      }

      //if this is a pit cell, move the water to the appropriate depression's water_vol.
      if(n == NO_FLOW){
        //Note that OCEAN cells have NO_FLOW so this will route water to the OCEAN
        //depression. This can be used to find out how much total run-off made it
        //into the ocean.

        if(arp.runoff(c)>0){
          deps[arp.label(c)].water_vol += arp.runoff(c)*arp.cell_area[y];
          //runoff is a depth of water and water_vol is a volume,
          //so multiply by the area of the pit cell to convert.
          assert(fp_ge(deps[arp.label(c)].water_vol,0) );
          if(deps[arp.label(c)].water_vol < 0)
            deps[arp.label(c)].water_vol = 0.0;
          arp.runoff(c) = 0; //Clean up as we go: this water is now stored in the DepressionHierarchy.
        }
      } else {                               //not a pit cell
        if(arp.runoff(c)>0){  //if there is water available, pass it downstream.
          //some infiltration happens as the water flows from cell to cell.
          //To do this, we need to first calculate the distance that the water travels.
          if(x == nx)      //same x means use e-w direction
            distance = arp.cellsize_e_w_metres[ny]/2.0;
          else if(y == ny) //this means use n-s direction
            distance = params.cellsize_n_s_metres/2.0;
          else    //it is diagonal, use pythagoras
            distance = (std::pow( (std::pow(arp.cellsize_e_w_metres[ny], 2.0) \
            + std::pow(params.cellsize_n_s_metres,2.0)), 0.5))/2.0;
          //divide each one by 2 because half of the infiltration will occur in the current cell,
          //and half in the neighbour cell receiving the water.

          //first, we do infiltration from the current cell, as the water moves
          //from the cell centre to the edge of the cell in the direction
          //of the neighbour.

          double my_infiltration = CalculateInfiltration(c,distance,static_cast<double>(arp.runoff(c)),arp);
          arp.wtd(c) += my_infiltration / static_cast<double>(arp.porosity(c));
          //add infiltration to the water table
          arp.runoff(c) -= my_infiltration;
          //and subtract the infiltration from the available surface water.
          assert(fp_le(arp.wtd(c),0));
          assert(fp_ge(arp.runoff(c),0));
          if(arp.wtd(c) > 0 )
            arp.wtd(c) = 0;
          if(arp.runoff(c) < 0)
            arp.runoff(c) = 0;

          //then, we do the same for the neighbour cell receiving the water,
          //from its edge that it receives along, to its centre.
          if(arp.runoff(c) > 0){
            //check again if there is water available since it's possible
            //the infiltration above used it up.
            my_infiltration = CalculateInfiltration(n,distance,static_cast<double>(arp.runoff(c)),arp);
            arp.wtd(n)    += my_infiltration / static_cast<double>(arp.porosity(n));
            arp.runoff(c) -= my_infiltration;


            assert(fp_le(arp.wtd(n),0));
            assert(fp_ge(arp.runoff(c),0));
            if(arp.wtd(n) > 0 )
              arp.wtd(n) = 0;
            if(arp.runoff(c) < 0)
              arp.runoff(c) = 0;
          }
        }

        //If we still have water, pass it downstream.
        if(arp.runoff(c)>0){
          //We use cell areas because the depth will change if we are moving
          //to a cell on a different latitude.
          arp.runoff(n) += arp.runoff(c)*arp.cell_area[y]/arp.cell_area[ny];
          //Add water to downstream neighbour.
          arp.runoff(c)  = 0;       //Clean up as we go: this water has been moved to the neighbouring cell.
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

  }
    RDLOG_TIME_USE<<"t FlowInDepressionHierarchy: Surface water = "<<timer.stop()<<" s";
}




///Here, we calculate the wtd_vol of all of the depressions.
///The reason for doing it here rather than in dephier is that I avoid
///having to make many adjustments to the wtd_vols during MoveWaterIntoPits.
///Wtd_vol represents the dep_vol plus any additional volume available to
///store water as groundwater.
///
///
///@param wtd          Water table depth. Values of 0 indicate saturation.
///                    Negative values indicate additional water can be added to the
///                    cell. Positive values indicate standing surface water.
///@param cell_area    The area of each cell, to calculate the total volume
///                    of water it can store.
///@param porosity     Porosity of the cell is needed to calculate the
///                    available below-ground water storage.
///@param final_label  Labels from GetDepressionHierarchy indicate which depression
///                    each cell belongs to. Final_label refers to the highest
///                    metadepression to which the cell belongs.
///@param deps         The DepressionHierarchy generated by GetDepressionHierarchy
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water storage space available (wtd_vol), including groundwater storage.
template<class elev_t, class wtd_t>
static void CalculateWtdVol(
  rd::Array2D<wtd_t>            &wtd,
  const std::vector<double>     &cell_area,
  const rd::Array2D<float>      &porosity,
  const rd::Array2D<dh_label_t> &final_label,
  DepressionHierarchy<elev_t>   &deps
){

  //We need to create yet another measure of volume, wtd_vol.
  //This will be the dep_vol plus the additional
  //volume allowed in the depression through storage as groundwater.
  //This is important because a depression may actually be able to store more
  //water than in its dep_vol. We may otherwise be overflowing a depression
  //when it is not actually supposed to overflow!
  //When we do overflow, we will also have to keep track of changes to the wtd
  //in the overflow depression, and associated changes in wtd_vol.
  //When a depression is completely saturated in groundwater, we will have wtd_vol == dep_vol.
  //In the ResetDH function, we have already set wtd_vol = dep_vol,
  //so here we only need to add any additional below-ground space.

  for(int y=0;y<wtd.height();y++)
  for(int x=0;x<wtd.width();x++){

    //cycle through the domain and add up all of the below-ground water storage space available
    auto clabel = final_label(x,y);
    //We should use final_label because a leaf depression will contain cells which are above its
    //outlet elevation. The groundwater space in these cells should not count towards the
    //wtd_vol of the leaf depression, but will be a part of a higher metadepression.
    //Final_label indicates the highest-level depression to which the cell belongs.

    if(clabel==OCEAN)  //We don't need to calculate ocean depression size.
      continue;

    assert(fp_le(wtd(x,y),0) );
    //because all wtds get set to max 0 when doing movewaterintopits -
    //surface water is now gathered in the pits and has not yet been distributed.
    if(wtd(x,y) > 0)
      wtd(x,y) = 0.0;

    //version for depth-variable porosity:
    //deps[clabel].wtd_only -= arp.cell_area[y] * arp.porosity(x,y) * arp.fdepth(x,y) * (exp(arp.wtd(x,y)/arp.fdepth(x,y)) - 1);
    deps[clabel].wtd_only -= static_cast<double>(cell_area[y]) * static_cast<double>(porosity(x,y)) * static_cast<double>(wtd(x,y));
    //this records the space that is only in the ground - not including above-ground space. This is needed to obtain the
    //appropriate total wtd_vols for metadepressions, below.
    //- because wtd is negative.
  }


  for(int d=0;d<(int)deps.size();d++){
    auto &dep = deps.at(d);
    if(dep.dep_label==OCEAN)
      continue;

    dep.wtd_vol = dep.dep_vol + dep.wtd_only;  //this records the total volume, both above and below ground, available to store water.


    if(dep.lchild!=NO_VALUE){  //if it has children, it is a metadepression and we need to add the groundwater space from the children.

      dep.wtd_vol      += deps.at(dep.lchild).wtd_only;     //store the total wtd_vols with all of your children included.
      dep.wtd_vol      += deps.at(dep.rchild).wtd_only;

      dep.wtd_only     += deps.at(dep.lchild).wtd_only;     //add this so that we know how much to add to the parent of this depression.
      dep.wtd_only     += deps.at(dep.rchild).wtd_only;
    }


    assert(fp_ge(dep.wtd_vol,0) );
    if(dep.wtd_vol < 0)
      dep.wtd_vol = 0.0;
    assert(fp_ge(dep.wtd_vol, dep.dep_vol) );  //wtd_vol cannot be smaller than dep_vol.
    if(dep.wtd_vol < dep.dep_vol)
      dep.wtd_vol = dep.dep_vol;
  }


  for(int d=0;d<(int)deps.size();d++){
    //Just checking that nothing went horribly wrong.
    auto &dep = deps.at(d);
    if(dep.dep_label==OCEAN)
      continue;
    if(dep.lchild!=NO_VALUE){
      assert(fp_ge(dep.dep_vol, deps.at(dep.lchild).dep_vol));
      assert(fp_ge(dep.dep_vol, deps.at(dep.rchild).dep_vol));
      assert(fp_ge(dep.wtd_vol, deps.at(dep.lchild).wtd_vol));
      assert(fp_ge(dep.wtd_vol, deps.at(dep.rchild).wtd_vol));
    }
  }
}



///We use this function to calculate any infiltration
///that occurs as water flows cell-by-cell downslope.
///This could happen both when the water is initially flowing down
///to the pit cells, and when it overflows from one depression to another.
///
///
///@param cell     The cell in which the infiltration will take place.
///@param distance The distance over which the water is travelling
///                (i.e. from one cell to the other)
///@param h_0      The amount of surface water available in the cell.
///
///@param arp      ArrayPack consisting of all of the various arrays used by this
///                and the other portions of the code. Here we use:
///                - arp.vert_ksat: The vertical hydraulic conductivity in the cell.
///                  This should differ from horizontal hydraulic conductivity,
///                  depending on local anisotropy. In the absence of better
///                  data, we recommend using soil types and a simple pedotransfer
///                  function to estimate this input.
///                - arp.slope: terrain slope, based on the topography data used.
///                - arp.wtd: The depth to water table. Here, we use it to check
///                  to see whether the cell is becoming fully saturated in
///                  groundwater during infiltration.
///                - arp.infiltration_array: This is an informational array,
///                  that records the total amount of infiltration that happens
///                  in each cell, should the user want this information.
///@return Calculates the amount of infiltration that will happen in this cell within
///        the time it takes the amount of  water to travel the given distance.
static inline double CalculateInfiltration(
  int         cell,
  double      distance,
  double      h_0,
  ArrayPack   &arp
){
  double mannings_n = 0.05; //TODO: is this the best value?
  double my_infiltration = 0.0;

 // double vert_ksat = std::max(static_cast<double>(arp.vert_ksat(cell)),0.000001);
   double vert_ksat = static_cast<double>(arp.vert_ksat(cell));
  //assigning a small minimum vertical ksat, since otherwise in cells with zero vert_ksat, the delta_t will be infinity.

  double slope = std::max(static_cast<double>(arp.slope(cell)), 0.000001);
  //assigning a small minimum slope, since otherwise in cells with zero slope, the delta_t will be infinity.


  //Manning's equation is rearranged to compute the amount of time it will take the given amount of water to cross the cell.
  double bracket = std::pow(h_0,(5.0/3.0)) - (vert_ksat * distance * (mannings_n/std::pow(slope,0.5)) * (5/3));
  double numerator = h_0 - std::pow(bracket,(3.0/5.0));
  double delta_t = numerator/vert_ksat;
  //delta_t is the amount of time it will take the given amount of
  //water to cross the given distance.
  if(vert_ksat == 0){
    return 0; //no infiltration can occur if there is zero vertical hydraulic conductivity.
  }

  //if it's a normal case, the water is not running out and the cell is not becoming saturated.
  //we calculate the infiltration using the delta_t from above, and ksat as a coefficient for infiltration amount.
  my_infiltration = vert_ksat*(delta_t);

  if(bracket < 0){           //all of the water is getting used up.
    my_infiltration =  h_0;  //We have to limit the infiltration to be equal to the available water, h0.
  }

  if( (arp.wtd(cell)*static_cast<double>(arp.porosity(cell)) )  > -my_infiltration){  //if the cell is getting saturated partway.
    my_infiltration = std::min(my_infiltration,static_cast<double>(-arp.wtd(cell))*static_cast<double>(arp.porosity(cell)) );
    //But, we must check for the case where both the cell becomes saturated
    //and the water gets used up.
    //sometimes the water may run out before the saturation occurred.
  }

  arp.infiltration_array(cell) += my_infiltration;

  return my_infiltration;
}




///At this point all the values of the water table `wtd`, which is not used by
///this function, are <=0. The excess water, which will eventually become
///standing surface water, is stored in the DepressionHierarchy itself in the
///leaf depressions.
///
///In this function we will perform a depth-first post-order traversal of the
///depression hierarchy starting with the OCEAN. When we reach the leaf
///depressions we check if `water_vol>wtd_vol`. If so, we try to overflow into
///the geographically proximal leaf depression indicated by our outlet. If there
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
///                           O(N) time. Only usable when infiltration is
///                           turned off.
///@param params              Various paramaters that are used here and in other
///                           portions of the code. Here, we don't use any of
///                           the parameters directly, but they are needed in
///                           OverflowInto which is called within this function.
///@param arp                 ArrayPack consisting of all of the various arrays used
///                           by this and the other portions of the code.
///                           The ArrayPack is not directly used in this function,
///                           but is used in MoveWaterinOverflow which is called
///                           by OverflowInto (which is in turn called by this)
///                           function, if infiltration is turned on.
///
///@return Modifies the depression hierarchy `deps` to indicate the amount of
///        water in each depression. This information can be used to add
///        standing surface water to the cells within a depression.
template<class elev_t>
static void MoveWaterInDepHier(
  int                                        current_depression,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  Parameters                                 &params,
  ArrayPack                                  &arp
){

  if(current_depression==NO_VALUE)
    return;

  auto &this_dep = deps.at(current_depression);

  //Catch depressions that link to the ocean through this one. These are special
  //cases because we will never spread water across the union of these
  //depressions and the current depressions: they only flow into the current
  //depression. We also process the ocean-linked depressions before the children
  //of this depression because the ocean-linked depressions can flow into the
  //children (but the children do not flow into the ocean-linked).
  for(const auto c: this_dep.ocean_linked)
    MoveWaterInDepHier(c, deps, jump_table,params,arp);

  //Visit child depressions. When these both overflow, then we spread water
  //across them by spreading water across their common metadepression
  MoveWaterInDepHier(this_dep.lchild, deps, jump_table, params, arp);
  MoveWaterInDepHier(this_dep.rchild, deps, jump_table, params, arp);

  //If the current depression is the ocean then at this point we've visited all
  //of its ocean-linked depressions (the ocean has no children). Since we do not
  //otherwise want to modify the ocean we quit here.
  if(current_depression==OCEAN)
    return;

  {
    const int lchild = this_dep.lchild;
    const int rchild = this_dep.rchild;

    //Only if both children are full should their water make its way to this
    //parent and, then, only if water hasn't been passed to the parent by
    //OverflowInto.
    if(lchild!=NO_VALUE  //If I am a leaf, then I can't get water from my children.
      && deps.at(lchild).water_vol>=deps.at(lchild).wtd_vol
      && deps.at(rchild).water_vol>=deps.at(rchild).wtd_vol
      && this_dep.water_vol == 0 //If water_vol>0 then children have already overflowed into parent
    ){
      this_dep.water_vol += deps.at(lchild).water_vol + deps.at(rchild).water_vol;
    }
    assert(fp_ge(this_dep.water_vol,0));
    if(this_dep.water_vol < 0)
      this_dep.water_vol = 0.0;
  }


  //Each depression has an associated dep_vol. This is the TOTAL volume of the
  //meta-depression including all of its children. This property answers the
  //question, "How much water can this meta-depression hold above ground?"
  //The wtd_vol then answers how much the depression can hold in TOTAL,
  //including below-ground.

  //Each depression also has an associated water volume called `water_vol`. This
  //stores the TOTAL volume of the water in the depression. That is, for each
  //depression this indicates how much water is in this depression and its
  //children. However, if a depression has sufficient volume to contain its
  //water, then its water_vol will not propagate up to its parent. In this way
  //we can distinguish between depressions whose water needs to be spread versus
  //metadepressions whose children might need water spread, but which will not
  //receive any spreading themselves.

  //We are overflowing the depression
  if(this_dep.water_vol>this_dep.wtd_vol) {



   if(current_depression==OCEAN)
    params.total_loss_to_ocean += this_dep.water_vol - this_dep.wtd_vol;

    //The neighbouring depression is not the ocean and this depression is
    //overflowing (therefore, all of its children are full)
    assert(this_dep.lchild==NO_VALUE || fp_eq(deps.at(this_dep.lchild).water_vol,deps.at(this_dep.lchild).wtd_vol));
    assert(this_dep.rchild==NO_VALUE || fp_eq(deps.at(this_dep.rchild).water_vol,deps.at(this_dep.rchild).wtd_vol));

    //OverflowInto will start at this depression and then, since there is
    //insufficient room, move to its neighbour via geolinks. If everything fills
    //up, the water will end up in this depression's parent. Afterwards, this
    //depression won't have extra water here any more. Any excess passed to the
    //parent will be dealt with it as we backtrack up the DH.

    OverflowInto(current_depression, -1, this_dep.parent, deps, jump_table, 0, params, arp);

    assert(
         fp_eq(this_dep.water_vol,0) || fp_le(this_dep.water_vol, this_dep.wtd_vol));

    assert((this_dep.lchild==NO_VALUE && this_dep.rchild==NO_VALUE)
      || ( this_dep.lchild!=NO_VALUE && this_dep.rchild!=NO_VALUE
        && fp_le(deps.at(this_dep.lchild).water_vol, this_dep.water_vol)
        && fp_le(deps.at(this_dep.rchild).water_vol, this_dep.water_vol)
         ));
  }

  //All overflowing depressions should by now have overflowed all the way down
  //to the ocean. We must now spread the water in the depressions by setting
  //appropriate values for wtd.
}




///This function is similar to MoveWaterIntoPits, it is also for routing
///water physically downslope and allows infiltration as it goes, where
///necessary.
///When it reaches the pit cell of a depression, all the remaining water is
///moved into the DepressionHierarchy data structure for rapid flood-spill-merge
///calculations.
///Because we are only moving water in one specific depression at a time from
///one specific inflow cell, I have used a while loop rather than a recursive
///call.
///This function will be called only when infiltration is turned on.
///
///@param extra_water   How much water, in total, is available to move downslope
///@param current_dep   The depression label we want to move water in
///@param previous_dep  The depression label the water came from
///@param deps          The DepressionHierarchy generated by GetDepressionHierarchy
///@param params        Various paramaters that are used here and in other portions
///                     of the code. Here we use:
///                     - params.cellsize_n_s_metres, to determine the distance
///                       travelled by water, if infiltration is on.
///@param arp           ArrayPack consisting of all of the various arrays used by this
///                     and the other portions of the code. Here we use:
///                     - arp.cell_area: Area of a cell, important for conserving
///                       water volume as it moves from cell to cell.
///                     - arp.cellsize_e_w_metres: cellsize in the east-west direction,
///                       used to calculated a distance travelled by flowing water
///                       when infiltration is turned on.
///                     - arp.flowdirs: the flow directions calculated in DH.
///                     - arp.label: the label of the leaf depression associated
///                       with each cell.
///                     - arp.final_label: the label of the top metadepression associated
///                       with each cell.
///                     - arp.topo: the topography, i.e. elevation above sea level
///                       of each cell.
///                     - arp.wtd: the depth to water table in each cell. If the water
///                       table is below the land surface (i.e. groundwater), then
///                       wtd will be negative. Positive values indicate the
///                       presence of surface water.
///
///@return  Modifies the wtd array to reflect infiltration which has occured.
///         Modifies the depression hierarchy `deps` to indicate the updated
///         amount of water contained in each depression.
template<class elev_t>
static void MoveWaterInOverflow(
  double                         extra_water,
  const dh_label_t               current_dep,
  const dh_label_t               previous_dep,
  DepressionHierarchy<elev_t>    &deps,
  Parameters                     &params,
  ArrayPack                      &arp

  ){
  auto &this_dep = deps.at(current_dep);
  int move_to_cell = NO_VALUE;  //the cell that will be the next to receive the extra water
  int previous_cell = deps.at(previous_dep).out_cell;  //This is the cell that the water is coming from.

  int my_label;

  int x,y;
  arp.topo.iToxy(previous_cell,x,y);
  double distance = 0;

  assert(extra_water > 0);

  int nx;  //the neighbour cells to check
  int ny;

  //The water must move downslope but we must also check that it moves
  //towards this_dep and doesn't flow back into the previous depression
  for(int n=1;n<=neighbours;n++){ //Check out our neighbours
    nx = x+dx[n];                 //Use offset to get neighbour x coordinate
    ny = y+dy[n];                 //Use offset to get neighbour y coordinate
    if(!arp.topo.inGrid(nx,ny))   //Is cell outside grid (too far North/South)?
      continue;                   //Yup: skip it.
    //we will choose this neighbour cell if it has the label of this_dep, so we are not going back into the
    //giving dep or any other adjacent deps. This is okay to do since we will only ever call
    //this function on a leaf depression: either when an ocean-linked depression overflows into
    //its leaf parent, or when we pass water to a geolinked neighbour.
    //We will also only choose this neighbour if its elevation is lower than any previously selected neighbour,
    //so that we can find the steepest downslope direction. We can't just use flowdirs here since we
    //need to ensure we choose a cell with the correct label.
    //TODO: Do we need to do anything specific for a case where there are multiple cells with the correct
    //label and identical elevations?

    if( ( this_dep.dep_label == arp.label(nx,ny)) && ( move_to_cell == NO_VALUE || arp.topo(nx,ny)<arp.topo(move_to_cell))){
      move_to_cell = arp.topo.xyToI(nx,ny);
    }

  }  //so now we know which is the correct starting cell.
  assert(move_to_cell != NO_VALUE);

  //we don't want to infiltrate in deps.at(previous_dep).out_cell if it does not share this depression's label,
  //because we have already done the overflow from that depression, so it could lead to an overfull depression if we do.
  //For simplicity, I am choosing not to infiltrate in the out_cell even if it shares this depression's label. TODO: fix this.

  while(extra_water > 0){
    assert(arp.label(move_to_cell) == this_dep.dep_label);

    int move_to_cell_x, move_to_cell_y;
    arp.topo.iToxy(move_to_cell,move_to_cell_x,move_to_cell_y);

    //Check the relative locations of the cells we are moving from and to, and then get the distance between them.
    //This distance is used in the call to CalculateInfiltration, below.
    if(x == move_to_cell_x)      //this means use e-w direction
      distance = arp.cellsize_e_w_metres[move_to_cell_y];  //distance travelled for infiltration
    else if(y == move_to_cell_y) //this means use n-s direction
      distance = params.cellsize_n_s_metres;
    else                         //diagonal, use pythagoras.
      distance = (std::pow( (std::pow(arp.cellsize_e_w_metres[move_to_cell_y], 2.0) + std::pow(params.cellsize_n_s_metres,2.0)), 0.5));

    //We must convert extra_water to a height when calling CalculateInfiltration.
    double my_infiltration = CalculateInfiltration(move_to_cell,distance,(extra_water/arp.cell_area[move_to_cell_y]),arp);

    if(my_infiltration > 0){
      //infiltration can't be greater than the available wtd space.
      assert(fp_ge(-(arp.wtd(move_to_cell) ), ( my_infiltration/ static_cast<double>(arp.porosity(move_to_cell)) ) ) );

      //If infiltration was greater than available wtd space by a tiny amount due to a precision error, fix that.
      if(-(arp.wtd(move_to_cell)) < (my_infiltration/ static_cast<double>(arp.porosity(move_to_cell)) )){
        my_infiltration = -arp.wtd(move_to_cell) * static_cast<double>(arp.porosity(move_to_cell));
        arp.wtd(move_to_cell) = 0.;
      }
      else{
        arp.wtd(move_to_cell) += my_infiltration / static_cast<double>(arp.porosity(move_to_cell));  //infiltrate into the cell
      }

      extra_water           -= my_infiltration *static_cast<double>(arp.cell_area[move_to_cell_y]) ;  //and that much of the water has been used up.

      assert(fp_le(arp.wtd(move_to_cell),0));  //we can't have a positive water table
      if(arp.wtd(move_to_cell)>=0)
        arp.wtd(move_to_cell) = 0;

      assert(fp_ge(extra_water,0));           //we can't have negative extra_water.

      //We have infiltrated into the cell and adjusted extra_water. Now, we need to adjust the wtd_vol of affected depressions.
      my_label = arp.label(move_to_cell);

      while(my_label != OCEAN ){      //adjust the wtd_vol of me and all my parents

        //only if it's below the outlet elevation do we modify wtd vol; cells higher on the slope never counted towards the wtd_vol to begin with.
        if(arp.topo(move_to_cell)<=deps.at(my_label).out_elev){
          deps.at(my_label).wtd_vol   -= my_infiltration*static_cast<double>(arp.cell_area[move_to_cell_y]);
        }

        assert(deps.at(my_label).wtd_vol - deps.at(my_label).dep_vol >= -FP_ERROR);
        if(deps.at(my_label).wtd_vol < deps.at(my_label).dep_vol)
          deps.at(my_label).wtd_vol = deps.at(my_label).dep_vol;

        if(deps.at(my_label).ocean_parent)  //Don't adjust the wtd_vols of any ocean parent depressions
          break;

        my_label = deps.at(my_label).parent;  //Go to the next parent in line to adjust its wtd_vol.
      }
      if(extra_water <= 0)
        break;
    }


    const auto ndir = arp.flowdirs(move_to_cell);
    if(ndir == NO_FLOW){   //we've reached a pit cell, so we can stop.
      this_dep.water_vol += extra_water;  //Add any remaining extra_water to this depression. If the depression is overfilled, it will still overflow in OverflowInto.
      extra_water = 0;    //No need to worry about infiltration in the pit cell: this will happen down in FillDepressions.
      break;
    }
    else{                 //we need to keep moving downslope.
      previous_cell = move_to_cell;
      arp.topo.iToxy(move_to_cell,x,y);
      //then we can use flowdirs to move the water all the way down to this_dep's pit cell.
      nx = x+dx[ndir];
      ny = y+dy[ndir];
      move_to_cell = arp.topo.xyToI(nx,ny);  //get the new move_to_cell
      assert(move_to_cell>=0);
    }
  }
}




///When water overflows from one depression into another, this function ensures
///that chained overflows and infiltration take place; that is, that water first
///fills the water table on the path between the two depressions and then moves
///from leaf depressions to metadepressions.
///
///A depression has three places it can put the water it's received.
///The depression will try to stash water in these three places sequentially.
///  1. It can store the water in itself
///  2. It can overflow into its neighbouring depression
///     (by following a geolink to that depression's leaf)
///  3. It can overflow into its parent
///
///Options (2) and (3) result in a recursive call. If there's enough water,
///eventually repeated calls to (3) will bring the function to the parent of the
///depression that originally called it (through its neighbour). At this point we
///stash the water in the parent and exit.
///
///Since we might end up calling the function again and again if there's a
///complex series of overflows, the `jump_table` argument holds the location of
///where the water ultimately ended up. Everything between the original node and
///this destination is then full which means that we only traverse each part of a
///filled hierarchy once.
///
///NOTE: This function is written recursively, but the jump table means it could
///follow a long string of depressions thereby hitting the stack/recursion limit.
///A better way of writing this would be to use an iterative stack technique, but
///this would also be slightly harder to describe, so we leave it for the future.
///
///@param root         The depression we're currently considering
///@param previous_dep The previous depression. When infiltration is turned on,
///                    we may need to know from which depression the water came.
///@param stop_node    When we reach this depression we dump all the excess water
///                    we're carrying into it. This depression is the parent of
///                    the depression that first called this function. Reaching
///                    it means that both the original metadepression and its
///                    neighbour are both full.
///@param deps         The DepressionHierarchy generated by
///                    GetDepressionHierarchy
///@param jump_table   A data structure that persists throughout the traversal
///                    and is used to skip from leaf nodes to the highest known
///                    meta-depression which still has unfilled volume. Ensures
///                    the traversal happens in O(N) time. Only usable when
///                    infiltration is set to off.
///@param extra_water  The amount of water left to distribute. We'll try to stash
///                    it in root. If we fail, we'll pass it to root's neighbour
///                    or, if the neighbour's full, to root's parent.
///@param params       Various paramaters that are used here and in other
///                    portions of the code. Here, we use:
///                    - params.infiltration_on, which tells us if infiltration
///                      is set to occur. If it is, we will sometimes call
///                      MoveWaterinOverflow when overflowing a depression.
///@param arp          ArrayPack consisting of all of the various arrays used
///                    by this and the other portions of the code.
///                    The ArrayPack is not directly used in this function,
///                    but is used in MoveWaterinOverflow which is called
///                    by this function, if infiltration is turned on.
///
///@return The depression where the water ultimately ended up. This is used to
///        update the jump table. The water_vol of various affected
///        depressions is changed within the body of the function.
template<class elev_t>
static dh_label_t OverflowInto(
  const dh_label_t                           root,
  const dh_label_t                           previous_dep,
  const dh_label_t                           stop_node,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,
  double                                     extra_water,
  Parameters                                 &params,
  ArrayPack                                  &arp
){

  auto &this_dep = deps.at(root);

  //If this depression is too full, then pick up its water. We'll find somewhere
  //to put it.
  if(this_dep.water_vol>this_dep.wtd_vol){
    //The marginal volume of this depression is larger than what it can hold, so
    //we determine the amount that overflows, the "extra water".
    extra_water += this_dep.water_vol - this_dep.wtd_vol;

    //Now that we've figured out how much excess water there is, we fill this
    //depression to its brim. Note that we don't use addition or subtraction
    //here to ensure floating-point equality.
    this_dep.water_vol = this_dep.wtd_vol;
  }


  //Determine if we've completed the loop. That'll happen either when we get to
  //the stop_node or the OCEAN. We can reach the OCEAN without reaching the
  //stop_node because the jump table can allow us to skip the stop_node.
  if(root==stop_node || root==OCEAN){
    //Otherwise, if we've reached the stop_node, that means we're at the
    //original node's parent. We give it our extra water.



    //in this case, then if infiltration is turned on, we should move water to the inlet of
    //this depression and let it flow downslope. Although this is only necessary if there
    //is groundwater space available AND if the water_vol is less than the wtd_vol.
    //If water_vol is equal to wtd_vol, then all of the groundwater space in the depression
    //will be filled later during FillDepressions in any case, so there is no
    //need to route water now.
    if(this_dep.wtd_vol > this_dep.dep_vol && this_dep.water_vol < this_dep.wtd_vol && params.infiltration_on == true){
      //okay, there is groundwater volume to fill, so we must move the water properly.

      //if this node is a normal parent, then I don't think it matters where in
      //the depression the water goes, and it can just get added to this_dep.water_vol.
      //But if the original depression is ocean_linked to this depression,
      //then the overflow is happening from a certain location and we need to route that water.
      //We also need to route the water if this depression was not the parent
      //of the previous depression (i.e. it was the odep). So, we check to see if
      //this was both the parent of the previous depression and if that was an ocean-linked
      //depression; or if this was not the parent of the previous depression.
      if( (  deps.at(previous_dep).parent       == this_dep.dep_label
          && deps.at(previous_dep).ocean_parent == true) \
          || deps.at(previous_dep).parent != this_dep.dep_label){


        MoveWaterInOverflow(extra_water,this_dep.dep_label,deps.at(previous_dep).dep_label,deps,params,arp);
      //  assert(this_dep.water_vol==0 || fp_le(this_dep.water_vol, this_dep.wtd_vol));
      }
    }
    else{
      this_dep.water_vol += extra_water;
      //no matter what, the water_vol needs to be updated.
      assert(fp_ge(this_dep.water_vol,0));
      if(this_dep.water_vol < 0){
        this_dep.water_vol = 0.0;
      }
    }

    return root;
  }

  //FIRST PLACE TO STASH WATER: IN THIS DEPRESSION

  if(this_dep.water_vol<this_dep.wtd_vol){                          //Can this depression hold any water?
    const double capacity = this_dep.wtd_vol - this_dep.water_vol;  //Yes. How much can it hold?

    if(extra_water >= capacity) {
      //It wasn't enough to hold all the water, so it will be completely filled
      //and there is no need to worry about flow routing.
      this_dep.water_vol = this_dep.wtd_vol;       //So we fill it all the way.
      extra_water       -= capacity;               //And have that much less extra water to worry about
    }
    else{  //extra_water < capacity, so we may have to worry about routing of
      //water and where it infiltrates. All of the extra_water will be added.

      //the depression has groundwater space, so we need to move water properly from the outlet.
      if(params.infiltration_on == true){
        if((this_dep.wtd_vol > this_dep.dep_vol) &&
          ((deps.at(previous_dep).parent == this_dep.dep_label && deps.at(previous_dep).ocean_parent == true)\
          || deps.at(previous_dep).parent != this_dep.dep_label)
         //but only if either I'm not the last dep's parent, or I am the parent but I'm an ocean-linking parent.
        ){
          MoveWaterInOverflow(extra_water,this_dep.dep_label,deps.at(previous_dep).dep_label,deps,params,arp);
        }
      }
      else
        this_dep.water_vol  = std::min(this_dep.water_vol + extra_water,this_dep.wtd_vol);
        //we should be using all of the extra water, but let's be careful about floating points.

      extra_water = 0;                        //No more extra water
    }
    assert(this_dep.water_vol==0 || fp_le(this_dep.water_vol, this_dep.wtd_vol));

  }

  if(fp_eq(extra_water,0))  {  //If there's no more extra water then call it quits
    return root;
  }

  //issue with the jump table: I can't figure out how to use MoveWaterInOverflow with it.
  //I require knowledge of the depression that sent the water to MoveWaterinOverflow, so that I know where to start moving the water downslope from.
  //When using the jump table, I know where the water ends up, but how do I know which side of the depression it flows in from?

  if(jump_table.count(root)!=0 && params.infiltration_on == false){  //we can only use the jump table if we dont' want water to infiltrate, or we won't know where infiltrating water is coming from.
    if(jump_table.at(root)==OCEAN)
      params.total_loss_to_ocean += extra_water;

    return jump_table[root] = OverflowInto(jump_table.at(root), -1, stop_node, deps, jump_table, extra_water, params, arp);
  }

  //Okay, so there's extra water and we can't fit it into this depression

  //SECOND PLACE TO STASH WATER: IN THIS DEPRESSION'S NEIGHBOUR
  //Maybe we can fit it into this depression's overflow depression!

  if(this_dep.odep != NO_VALUE){  //Does the depression even have such a neighbour?
    auto &odep = deps.at(this_dep.odep);
    if(odep.water_vol<odep.wtd_vol){    //Can overflow depression hold more water?
      //Yes. Move the water geographically into that depression's leaf.
      if(this_dep.geolink==OCEAN)
        params.total_loss_to_ocean += extra_water;
      return jump_table[root] = OverflowInto(this_dep.geolink, this_dep.dep_label, stop_node, deps, jump_table, extra_water, params, arp);
    } else if(odep.water_vol >= odep.wtd_vol) {
      //Neighbour is overfull. Since our next stop is the parent, let's take the
      //neighbour's water with us when we go.
      extra_water += odep.water_vol-odep.wtd_vol;
      odep.water_vol = odep.wtd_vol;
    }
  }

  //Okay, so the extra water didn't fit into this depression or its overflow
  //depression. That means we pass it to this depression's parent.Since we're
  //going to add water to that parent, we need it to know about all the water
  //beneath it, so we add our water and our overflow neighbour's water before
  //climbing up to the parent. Since it's possible our overflow neighbour has
  //already done this, we only add the water we and our neighbour own if the
  //parent currently has no water in it. If we're ocean-linked to our parent we
  //don't want to do this since we are not contained within our parent.
  auto &pdep = deps.at(this_dep.parent);
  if(pdep.water_vol==0 && !this_dep.ocean_parent){
    pdep.water_vol += this_dep.water_vol;
    if(this_dep.odep!=NO_VALUE){
      pdep.water_vol += deps.at(this_dep.odep).water_vol;
    }
  }




  //THIRD PLACE TO STASH WATER: IN THIS DEPRESSION'S PARENT
     if(this_dep.parent==OCEAN)
      params.total_loss_to_ocean += extra_water;
  return jump_table[root] = OverflowInto(this_dep.parent, this_dep.dep_label, stop_node, deps, jump_table, extra_water, params, arp);
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
  int   leaf_label = NO_VALUE;
  //The metadepression containing all of the children. This metadepression is
  //guaranteed to be large enough to hold all of the water of its children plus
  //whatever exists only in the metadepression itself. We use this to determine
  //the water and depression volumes.
  int   top_label = NO_VALUE;
  //Here we keep track of which depressions are contained within the
  //metadepression. This allows us to limit the spreading function to cells
  //within the metadepression.
  std::unordered_set<dh_label_t> my_labels;
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
///@param arp          ArrayPack consisting of all of the various arrays used
///                    by this and the other portions of the code.
///                    The ArrayPack is not directly used in this function,
///                    but is used in FillDepressions and FillAFullDepression
///                    which are called by this function.
///
///@return Information about the subtree: its leaf node, depressions it
///        contains, and its root node.
template<class elev_t>
static SubtreeDepressionInfo FindDepressionsToFill(
  const int                         current_depression,  //Depression we are currently in
  DepressionHierarchy<elev_t>       &deps,               //Depression hierarchy
  ArrayPack                         &arp
){
  //Stop when we reach one level below the leaves
  if(current_depression==NO_VALUE)
    return SubtreeDepressionInfo();

  const auto& this_dep = deps.at(current_depression);

  //We start by visiting all of the ocean-linked depressions. They don't need to
  //pass us anything because their water has already been transferred to this
  //metadepression tree by MoveWaterInDepHier(). Similarly, it doesn't matter
  //what their leaf labels are since we will never spread water into them.
  for(const auto c: this_dep.ocean_linked)
    FindDepressionsToFill(c, deps, arp);

  //At this point we've visited all of the ocean-linked depressions. Since all
  //depressions link to the ocean and the ocean has no children, this means we
  //have visited all the depressions and spread their water. Since we don't wish
  //to modify the ocean, we are done.
  if(current_depression==OCEAN)
    return SubtreeDepressionInfo();

  //We visit both of the children. We need to keep track of info from these
  //because we may spread water across them.
  SubtreeDepressionInfo left_info  = FindDepressionsToFill(this_dep.lchild, deps, arp);
  SubtreeDepressionInfo right_info = FindDepressionsToFill(this_dep.rchild, deps, arp );

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
  assert(fp_le(this_dep.water_vol, this_dep.wtd_vol));

  //If the depression holds less water than its volume, we can fill it now.
  //Alternatively, if its parent is ocean-linked then we need to fill the
  //depression now because the excess will have already overflowed through the
  //oceanlink and we won't get another chance (oceanlinked depressions may be at
  //the bottom of a cliff with respect to the focal depression). Alternatively,
  //if the depression's parent contains no water then we know the depression is
  //large enough to contain its water.
  if(this_dep.water_vol<this_dep.wtd_vol || this_dep.ocean_parent || \
  (this_dep.water_vol==this_dep.wtd_vol && deps.at(this_dep.parent).water_vol==0)){
    assert(fp_le(this_dep.water_vol, this_dep.wtd_vol));


    //If the depression is exactly full, then we already know the level
    //that the water should be: the outlet of the depression. Due to large
    //depression sizes, the methods in FillDepressions can have difficulty
    //with precision in cases where the depression is exactly full.
    //Therefore, we use FillAFullDepression to check the exact cells contained
    //in this depression and fill them to the level that we already know.
    if(this_dep.water_vol == this_dep.wtd_vol){
      FillAFullDepression(combined,deps,arp);
    }

    //If both of a depression's children have already spread their water,
    //we do not want to attempt to do so again in an empty parent depression.
    //We check to see if both children have finished spreading water.
    else{
      FillDepressions(deps.at(combined.leaf_label).pit_cell,
                      deps.at(combined.top_label).out_cell,
                      deps.at(combined.top_label).out_elev,
                      combined.my_labels, this_dep.water_vol, deps,arp);
    }
    //At this point there should be no more water all the way up the tree until
    //we pass through an ocean link, so we pass this up as a kind of null value.
    return SubtreeDepressionInfo();
  } else {
    return combined;
  }
}


///This function adjusts the water table to reflect standing surface water that
///has pooled at the bottom of depressions, in the special case when
///the depression is exactly full, i.e. the water volume equals the
///wtd_vol of the depression. This means that we already know the water level:
///the elevation of the outlet cell. We need only find which cells are
///contained within this depression to update the water level within these
///cells.
///
///We know the ID of the depression to fill, as well as the labels of
///any smaller depressions contained within it.
///
///All leaf depressions store within their DH metadata a list of which cells
///they contain. By cycling though the cells contained wihtin any
///affected depressions, we can quickly locate the cells affected
///by this water level update.
///
///@param stdi     Leaf node, metadepression, depressions between. Determines
///                the extent of the flooding.
///@param deps     The DepressionHierarchy generated by GetDepressionHierarchy
///@param topo     Topography used to generate the DepressionHierarchy
///@param label    Labels from GetDepressionHierarchy indicate which depression
///                each cell belongs to.
///@param arp      ArrayPack consisting of all of the various arrays used by this
///                and the other portions of the code. Here we use:
///                - arp.topo: we need to know the elevations of cells within
///                  the affected depressions.
///                - arp.wtd: The depth to water table. Here, we update this
///                  to reflect the new water level within the depression
///                  being filled.
///
///@return         N/A
template<class elev_t>
static void FillAFullDepression(
    SubtreeDepressionInfo             &stdi,
    DepressionHierarchy<elev_t>       &deps,
    ArrayPack                         &arp
){
  //scan through all of the labels which are a part of this metadepression
  for(const auto label: stdi.my_labels){
    //for each one, check through all of the cells contained within the
    //depression with that label (only leaf depressions store this information)
    for(const auto cell: deps.at(label).my_cells){
      //if the cell is lower than the outlet elevation, it is affected
      //by the depression fill on this level, so we update the wtd in this cell.
      if(arp.topo(cell) > deps.at(stdi.top_label).out_elev)
        break;
      else
        arp.wtd(cell) = deps.at(stdi.top_label).out_elev - static_cast<double>(arp.topo(cell));
    }
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
///@param pit_cell   Cell to start trying to fill from. This is the pit cell
///                  of one of the leaf depressions within this depression.
///@param out_cell   The outlet cell of this depression.
///@param out_elev   The elevation of the outlet cell of this depression.
///@param dep_labels Labels contained within the metadepression we are trying to
///                  fill
///@param water_vol How much water needs to be spread throughout the depression
///@param deps      The DepressionHierarchy generated by GetDepressionHierarchy
///@param arp       ArrayPack consisting of all of the various arrays used by this
///                 and the other portions of the code. Here we use:
///                 - arp.topo: we need to know the elevations of cells within
///                   the affected depressions.
///                 - arp.label: we check labels to see if a cell belongs to
///                   the affected depressions.
///                 - arp.cell_area: we need to compute the volume within
///                   a given portion of the depression, requiring the areas
///                   of the cells it contains.
///                 - arp.porosity: required because we will fill in the
///                   water table when below the land surface as we go.
///                 - arp.wtd: The depth to water table. Here, we update this
///                   to reflect the new water level within the depression
///                   being filled.
///@return         N/A
template<class elev_t>
static void FillDepressions(
  const int                            pit_cell,
  const int                            out_cell,
  const double                         out_elev,
  const std::unordered_set<dh_label_t> &dep_labels,
  double                               water_vol,
  DepressionHierarchy<elev_t>          &/*deps*/,
  ArrayPack                            &arp
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
  std::unordered_set<flat_c_idx> visited(2048);

  //Priority queue that sorts cells by lowest elevation first. If two cells are
  //of equal elevation the one added most recently is popped first. The ordering
  //of the cells processed by this priority queue need not match the ordering of
  //the cells processed by the depression hierarchy.
  GridCellZk_high_pq<elev_t> flood_q;

  //Cell from which we can begin flooding the meta-depression. Which one we
  //choose is arbitrary, since we will fill all of the leaf depressions and
  //intermediate metadepressions until and including when we reach the
  //depression identified by stdi.top_label
  assert(pit_cell>=0);

  //We start flooding at the pit cell of the depression and work our way
  //upwards
  flood_q.emplace(
    pit_cell % arp.topo.width(),
    pit_cell / arp.topo.width(),
    arp.topo(pit_cell)
  );

  visited.emplace(pit_cell);


  //Cells whose wtd will be affected as we spread water around
  std::vector<flat_c_idx> cells_affected;

  //Stores the sum of the elevations of all of the cells in cells_affected. Used
  //for calculating the volume we've seen so far. (See explanation above or in
  //dephier.hpp TODO)
  double current_elevation = 0;   //elevation of the cell we are now examining
  double previous_elevation = 0;  //elevation of the cell examined just before this
  double area_times_elevation_total = 0; //used to compute water level
  double current_volume = 0;      //volume in the depression based on the cells visited so far
  double current_area = 0;        //total area of the cells visited so far

  while(!flood_q.empty()){
    const auto c = flood_q.top();
    flood_q.pop();
    current_elevation = static_cast<double>(arp.topo(c.x,c.y));

    //We keep track of the current volume of the depression by noting the elevation
    //of the cell we are currently examining, that of the previous cell, and the
    //total area of the depression examined so far.

    //Note that the current cell's above ground volume and wtd do not
    //contribute at all. It is as though there is a virtual water line
    //coincident with the edge of the current cell.
    //No water infiltrates into this cell or is stored above it -
    //only cells previously visited are considered when doing volume
    //calculations.

    //Current volume of this subset of the metadepression. Since we might climb
    //over a saddle point, this value can occasionally be negative. It will be
    //positive by the time we need to spread the water around.

    current_volume += (current_elevation-previous_elevation)*current_area;

    assert(fp_ge(water_vol,0) );
    if(water_vol < 0)
      water_vol = 0;

    //All the cells within this depression should have water table depths less
    //than or equal to zero because we have moved all of their water down slope
    //into the pit cell. Since we may already have filled other depressions
    //their cells are allowed to have wtd>0. Thus, we raise a warning if we are
    //looking at a cell in this unfilled depression with wtd>0.
    if(dep_labels.count(arp.label(c.x,c.y))==1 && arp.wtd(c.x,c.y)>FP_ERROR)
      throw std::runtime_error("A cell was discovered in an unfilled depression with wtd>0!");

    //There are two possibilities:
    //1. The virtual water level exceeds slightly the height of the cell.
    //The cell's water table then fills up as much as it can.
    //   The water surface is then level with the height of the cell.
    //
    //2. The water surface is below the height of the cell because there is
    //sufficient topographic volume to hold all the water.
    //   In this case, the cell's water table is left unaffected.


//depth-variable porosity version:
//    if(water_vol - (current_volume - (arp.cell_area[c.y] * arp.porosity(c.x,c.y) * arp.fdepth(c.x,c.y) * (exp(arp.wtd(c.x,c.y)/arp.fdepth(c.x,c.y)) - 1) ) ) <= FP_ERROR){


    if( (((water_vol - (current_volume - (static_cast<double>(arp.cell_area[c.y]) * static_cast<double>(arp.porosity(c.x,c.y)) * \
      static_cast<double>(arp.wtd(c.x,c.y)) ) ) )<=FP_ERROR) && dep_labels.count(arp.label(c.x,c.y))==1) \
      || water_vol-current_volume < FP_ERROR ) {
      //we are checking for two similar but different cases here: either the current volume is large enough to accommodate
      //all of the water, without needing to check for groundwater space in the current cell; or
      //there is enough space if we account for groundwater space in the current cell AND the current
      //cell is labelled as a part of this depression.
      //The distinction is necessary because the outlet cell of the depression does not always
      //share a label with this depression. When that happens, we do not want to allow infiltration
      //in that outlet cell since it contributes to the wtd_vol of a different depression.

      //The current scope of the depression plus the water storage capacity of
      //this cell is sufficient to store all of the water. We'll stop adding
      //cells and fill things now.

      auto water_level = DetermineWaterLevel(arp.topo,arp.porosity,arp.cell_area,arp.wtd, water_vol, c.x,  c.y, current_volume, current_area, area_times_elevation_total);

      //Water level must be higher than (or equal to) the previous cell we looked at, but lower than (or equal to) the current cell
      assert(cells_affected.size()==0 || fp_le(arp.topo(cells_affected.back()),water_level) );
      assert(fp_ge(arp.topo(c.x,c.y),water_level));

      BackfillDepression(water_level, arp.topo, arp.wtd, cells_affected);

      //We've spread the water, so we're done
      return;

    }

    //We haven't found enough volume for the water yet.

    //The outlet cell never impacts the amount of water the depression can hold
    //and its water table is adjusted in the if-clause above, if appropriate.
    if(static_cast<int32_t>(arp.topo.xyToI(c.x,c.y)) != out_cell ){

      //Add this cell to those affected so that its volume is available for
      //filling.
      cells_affected.emplace_back(arp.topo.xyToI(c.x,c.y));

      //Fill in cells' water tables as we go
      assert(fp_le(arp.wtd(c.x,c.y),0) );
      if(arp.wtd(c.x,c.y) > 0)
        arp.wtd(c.x,c.y) = 0;


   //  depth-variable porosity version
   //   water_vol += arp.cell_area[c.y] * arp.porosity(c.x,c.y) * arp.fdepth(c.x,c.y) * (exp(arp.wtd(c.x,c.y)/arp.fdepth(c.x,c.y)) - 1);
      water_vol +=  static_cast<double>(arp.cell_area[c.y]) * static_cast<double>(arp.porosity(c.x,c.y)) * static_cast<double>(arp.wtd(c.x,c.y)) ;
      //We use += because wtd is less than or equal to zero
      arp.wtd(c.x,c.y)    = 0;      //Now we are sure that wtd is 0, since we've just filled it

      //Add the current cell's information to the running total
      current_area += static_cast<double>(arp.cell_area[c.y]);
      //adding to the area after the volume because when there is only 1 cell,
      //the answer for volume should be 0, etc.
      //We don't want to include the area of the target cell.
      area_times_elevation_total += arp.topo(c.x,c.y)*arp.cell_area[c.y];
    }

    //Add focal cell's neighbours to expand the search for more volume.
    for(int n=1;n<=neighbours;n++){
      const int nx = c.x + dx[n];            //Get neighbour's x-coordinate using an offset and wrapping
      //TODO ModFloor(x+dx[n],topo.width());
      const int ny = c.y + dy[n];            //Get neighbour's y-coordinate using an offset
      if(!arp.topo.inGrid(nx,ny))            //Is this cell in the grid?
        continue;                            //Nope: out of bounds.
      const int ni = arp.topo.xyToI(nx,ny);  //Get neighbour's flat index

      //Don't add cells which are not part of the depression unless the cell in
      //question is the outlet cell.
      if(dep_labels.count(arp.label(ni))==0 && ni!=out_cell)
        continue;

      //Don't add cells which are higher than the outlet elevation:
      //if these cells are being reached, it could be time to add
      //the outlet cell to the priority queue.
      if(arp.topo(nx,ny) > out_elev)
        continue;

      if(visited.count(ni)==0){
        //add any neighbour cells which have not yet been checked, to the priority queue.
        flood_q.emplace(nx,ny,arp.topo(nx,ny));
        visited.emplace(ni);
      }
    }

    //The queue is empty, so we add the outlet cell. We may have already visited
    //the cell while moving from one side of a depression to the other, but
    //since it doesn't affect depression volume it's safe to visit it twice.
    //NOTE: if our logic is wrong somewhere adding the cell like this could
    //result in an infinite loop.
    if(flood_q.empty()){
      int x,y;
      arp.topo.iToxy(out_cell , x, y);
      flood_q.emplace(x, y, arp.topo(out_cell));
      visited.emplace(out_cell);
    }

    previous_elevation = static_cast<double>(arp.topo(c.x,c.y));
  }

  //Since we're in this function we are supposed to be guaranteed to be able to
  //fill our depression, since we only enter this function once that is true.
  //Therefore, if we've reached this point, something has gone horribly wrong
  //somewhere. :-(

  assert(!flood_q.empty());
  throw std::runtime_error("PQ loop exited without filling a depression!");
}


///Once a set of cells sufficient to hold the water in a depression has been
///found, we use this function to spread the water across those cells.
///@param water_level    Elevation of the water's surface
///@param topo           Elevations of cells
///@param wtd            Water table depth
///@param cells_affected Queue of cells to add water to
template<class elev_t, class wtd_t>
void BackfillDepression(
  double                  water_level,
  const Array2D<elev_t>   &topo,
  Array2D<wtd_t>          &wtd,
  std::vector<flat_c_idx> &cells_affected
){
  for(const auto c: cells_affected){
    //While searching for cells to fill, we simulated infiltration. The result
    //is that surface water should now be covering cells whose water table is
    //flush with the surface.
    assert(fp_ge(wtd(c),0));
    if(wtd(c) < 0)
      wtd(c) = 0;

    //Water level should be equal to or greater than the elevation of each cell
    //to be filled
    assert(fp_ge(water_level,topo(c)) );
    if(water_level < topo(c))
      water_level = topo(c);

    //Raise the cell's water table so that it is flush with the water's surface
    //(Note that the elevation of the surface is topo+wtd)
    wtd(c) = water_level - topo(c);
    //The water table must have been zero when we added the water, but there's a
    //chance floating point effects will make it negative. Therefore, we make
    //sure that it's >=0. This is most likely to be necessary when considering
    //saddle cells within a metadepression.
    wtd(c) = std::max(static_cast<wtd_t>(0), wtd(c));
    assert(wtd(c)>=0);
  }
}



///The Lake-Level Equation.
///@param sill_wtd               Water table depth of the sill cell
///@param water_vol              Volume of water to spread across the depression
///@param sill_elevation         Elevation of the cell that dams the water into
///                              the depression.
///@param cells_to_spread_across How many cells the water is to be spread across
///
///@return The elevation of the water level
template<class elev_t, class wtd_t>
double DetermineWaterLevel(
  const Array2D<elev_t>     &topo,
  const Array2D<float>      &porosity,
  const std::vector<double> &cell_area,
  Array2D<wtd_t>            &wtd,
  double       water_vol,
  const int    cx,
  const int    cy,
  const double current_volume,
  const double current_area,
  const double area_times_elevation_total
){
  double water_level;
  if(current_volume<water_vol){
    //The volume of water exceeds what we can hold above ground, so we will
    //stash as much as we can in this cell's water table. This is okay
    //because the above ground volume plus this cell's water table IS enough
    //volume (per the if-clause above).

    //Fill in as much of this cell's water table as we can
    const double fill_amount = water_vol - current_volume;
    assert(fp_ge(fill_amount,0));

    //depth-variable porosity version:
    //    arp.wtd(c.x,c.y) = arp.fdepth(c.x,c.y) * log(exp(arp.wtd(c.x,c.y)/arp.fdepth(c.x,c.y)) +
    //      fill_amount/(arp.cell_area[c.y] * arp.porosity(c.x,c.y) * arp.fdepth(c.x,c.y)));

    wtd(cx,cy) += fill_amount/cell_area[cy]/porosity(cx,cy);

    water_level = topo(cx,cy);


  } else if (current_volume==water_vol) {
    //The volume of water is exactly equal to the above ground volume
    //so we set the water level equal to this cell's elevation
    water_level = topo(cx,cy);
  } else {  //The water volume is less than this cell's elevation,
    //so we calculate what the water level should be.
    //Volume of water = sum(current_area*(water_level -
    //cell_elevation))
    //                = sum(current_area * water_level) -
    //sum(current_area * cell_elevation)
    //                = water_level * sum(current_area) -
    //sum(current_area * cell_elevation)
    //rearranging this gives us:

    water_level = (water_vol / current_area) + (area_times_elevation_total / current_area);
  }
  return water_level;
}

///Reset water volumes in the Depression Hierarchy to zero
///
///@param deps Depression Hierarchy to reset
template<class elev_t>
void ResetDH(DepressionHierarchy<elev_t> &deps){
  for(auto &dep: deps){
    dep.water_vol = 0;
    dep.wtd_only = 0;  //used later in wtd_vol calculations
  }
}

}
