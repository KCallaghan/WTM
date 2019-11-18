#include "transient_groundwater.hpp"
#include "fill_spill_merge.hpp"
#include "evaporation.hpp"

#include "../common/netcdf.hpp"
#include "ArrayPack.hpp"
#include "parameters.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include <fstream>
using namespace std;
 


namespace rd = richdem;
namespace dh = richdem::dephier;

const double UNDEF  = -1.0e7;


int main(int argc, char **argv){

  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <Configuration File>"<<std::endl;
    return -1;
  }

  ArrayPack arp;
  std::cerr<<"Argv"<<argv<<std::endl;
  Parameters params(argv[1]);


//Initialise all of the values for an transient-style run:
//Open data files, set the changing cellsize arrays, do any needed unit conversions. 

  ofstream textfile;
  textfile.open ("text_coupled_transient.txt", std::ios_base::app);



  textfile<<"Initialise transient"<<std::endl;



 //load in the data files: ksat, mask, e-folding depth, precipitation, topography, and evaporation files. 


//TODO: Load in 'end' values for each parameter, for transient runs
  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");   //TODO: check units of ksat. I think it is m/s. 

  params.ncells_x = arp.ksat.width();  //width and height in number of cells in the array
  params.ncells_y = arp.ksat.height();

  arp.land_mask     = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); //A binary mask that is 1 where there is land and 0 in the ocean
  arp.fdepth        = LoadData<float>(params.surfdatadir + params.region + "_fdepth_fully_calculated.nc", "value");  //fslope = 100/(1+150*slope), f>2.5 m. Note this is specific to a 30 arcsecond grid! Other grid resolutions should use different constants. 
  //TODO: Note the previsor that f>2.5 m! I don't think I have done this in the past, check! 
  arp.precip        = LoadData<float>(params.surfdatadir + params.region + "_precip.nc",   "value");
  arp.temp          = LoadData<float>(params.surfdatadir + params.region + "_temp.nc",   "value");
  arp.topo          = LoadData<float>(params.surfdatadir + params.region + "_topo.nc",   "value");
  arp.starting_evap = LoadData<float>(params.surfdatadir + params.region + "_evap.nc",   "value");
  arp.relhum        = LoadData<float>(params.surfdatadir + params.region + "_relhum.nc",   "value");

  //Initialise arrays that start out empty:

  arp.wtd    = rd::Array2D<float>(arp.topo,0.0f);  //TODO: load in an actual wtd array. This should come from previous model runs and not start out as 0.
  arp.rech   = rd::Array2D<float>(arp.topo,0);  //We have to calculate recharge within the code since we will change the way evaporation works depending on whether or not there is surface water present. 
  arp.head   = rd::Array2D<float>(arp.topo,0);        //Just to initialise these - we'll add the appropriate values later. 
  arp.kcell  = rd::Array2D<float>(arp.topo,0);
  arp.evap   = rd::Array2D<float>(arp.topo,0);  
  arp.e_sat  = rd::Array2D<float>(arp.topo,0);  
  arp.e_a    = rd::Array2D<float>(arp.topo,0);  

  arp.wtd_change_total = rd::Array2D<float>(arp.topo,0.0f);   //This array is used to store the values of how much the water table will change in one iteration, then adding it to wtd gets the new wtd.
  arp.stability_time_seconds = rd::Array2D<float>(arp.topo,0);  //This is used to check whether the calculation is stable for the selected time step. TODO: Should I just remove this array if we're not going to do adaptive time-stepping?


  textfile<<"loaded all"<<std::endl;

  //compute changing cell size and distances between cells as these change with latitude:

  const float earth_radius = 6371000.; //metres

  params.cellsize_n_s_metres = (earth_radius*(M_PI/180.))/params.cells_per_degree; //distance between lines of latitude is a constant. 

//initialise some arrays
  arp.latitude_radians.resize      (params.ncells_y);   // the latitude of each row of cells
  arp.cellsize_e_w_metres.resize  (params.ncells_y);    //size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cellsize_e_w_metres_N.resize  (params.ncells_y);  //size of a cell in the east-west direction at the northern edge of the cell (metres)
  arp.cellsize_e_w_metres_S.resize  (params.ncells_y);  //size of a cell in the east-west direction at the southern edge of the cell (metres)
  arp.cell_area.resize (params.ncells_y);               //cell area (metres squared)


  for(unsigned int j=0;j<arp.latitude_radians.size();j++){
    //latitude at the centre of a cell:
    arp.latitude_radians[j]       = (float(j)/params.cells_per_degree+params.southern_edge)*(M_PI/180.); //southern edge of the domain in degrees, plus the number of cells up from this location/the number of cells per degree, converted to radians.

    //cells_per_degree = 120, there are this many 30 arc-second pieces in one degree. 
    //j/cells_per_degree gives the number  of degrees up from the southern edge, add southern_edge since the southern edge may not be at 0 latitude. *pi/180 to convert to radians. 
    //latitude_radians is now the latitude in radians. 

    //latitude at the southern edge of a cell (subtract half a cell):
    double latitude_radians_S = ((float(j) - 0.5)/params.cells_per_degree+params.southern_edge)*(M_PI/180.); 
    //latitude at the northern edge of a cell (add half a cell):
    double latitude_radians_N = ((float(j) + 0.5)/params.cells_per_degree+params.southern_edge)*(M_PI/180.); 

    //distance between lines of longitude varies with latitude. This is the distance at the centre of a cell for a given latitude:
    arp.cellsize_e_w_metres[j] = earth_radius*std::cos(arp.latitude_radians[j])*(M_PI/180.)/params.cells_per_degree;

    //distance at the northern edge of the cell for the given latitude:
    arp.cellsize_e_w_metres_N[j] = earth_radius*std::cos(latitude_radians_N)*(M_PI/180.)/params.cells_per_degree;
    //distance at the southern edge of the cell for the given latitude:
    arp.cellsize_e_w_metres_S[j] = earth_radius*std::cos(latitude_radians_S)*(M_PI/180.)/params.cells_per_degree;

    //cell area computed as a trapezoid, using unchanging north-south distance, and east-west distances at the northern and southern edges of the cell:
    arp.cell_area[j] = params.cellsize_n_s_metres* (arp.cellsize_e_w_metres_N[j] + arp.cellsize_e_w_metres_S[j])/2;


    delete[] arp.cellsize_e_w_metres_S;
    delete[] arp.cellsize_e_w_metres_N;
    delete[] arp.latitude_radians;
//TODO: see which, if any, arrays can be cleared after use to free up memory. How to do this?


  }

  textfile<<"computed distances, areas, and latitudes"<<std::endl;

  //Change undefined cells to 0
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo(i)<=UNDEF){
      arp.topo(i) = 0;
    }
    if(arp.ksat(i)<=0){
      arp.ksat(i) = 1E-10; //ksat is always supposed to be greater than 0. This is a small value (impermeable material)

      //TODO: I have now exported ksat to be only >= 1E-10 so can remove this
      //TODO: ksat is in m/s, do I need to convert units for the time step?
    }
    //Converting to appropriate time step
    arp.rech(i) *= 100;           //TODO: This should be converted to whatever the timestep is, not necessarily monthly. It should check the time step from the params.
    
    //do we need to do a unit conversion (m -> mm) for rech? For equilibrium we set values >10000 to 0, do we need to do this?

  } 


  arp.check();

//get the starting recharge using precip and evap inputs:


//TODO: Note these units are m/yr. convert as needed. 
//TODO: Use actual evap layers to get the appropriate precip, this is a placeholder. 
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.evap(i) = arp.starting_evap(i);
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);
  }


//Wtd is 0 in the ocean:
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo.isNoData(i) || arp.topo(i)==dh::OCEAN_LEVEL || arp.land_mask(i) == 0.0f){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      arp.wtd  (i) = 0;
    }
  }

 //We have finished initialising all of the arrays and data that will be used for groundwater. 


//initialise the things we need for surface water:

  rd::Array2D<dh::dh_label_t> label   (params.ncells_x, params.ncells_y, dh::NO_DEP ); //No cells are part of a depression
  rd::Array2D<dh::dh_label_t> final_label   (params.ncells_x, params.ncells_y, dh::NO_DEP ); //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(params.ncells_x, params.ncells_y, rd::NO_FLOW); //No cells flow anywhere


 //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(arp.topo.isNoData(i) || arp.topo(i)==dh::OCEAN_LEVEL || arp.land_mask(i) == 0.0f){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      label(i) = dh::OCEAN;
      final_label(i) = dh::OCEAN;
      arp.wtd  (i) = 0;
    }
  }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(arp.topo, label, final_label, flowdirs);



int cycles_done = 0;

 float max_total = 0.0;
 float min_total = 0.0;

while(true){

  textfile<<"Cycles done: "<<cycles_done<<std::endl;
  if(cycles_done == 500000)  //How many times to repeat the process
    break;

  if((cycles_done % 10) == 0){
    textfile<<"saving partway result"<<std::endl;  
    SaveAsNetCDF(arp.wtd,"N_America_transient-surface-coupled-full-res.nc","value");

  }

 //Run the transient groundwater code to move water


//TODO: Update data arrays every x number of years so they change from starting to end values through time
  arp.wtd = transient(params,arp);



//Move surface water.
  dh::FillSpillMerge(arp.topo, label, final_label, flowdirs, deps, arp.wtd);


  evaporation_update(params,arp);


  cycles_done += 1;
}

//We are finished, save the result.   
  textfile<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"N_America_transient-surface-coupled-full-res.nc","value");

  textfile.close();


  return 0;

}
