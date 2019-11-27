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

  if(argc!=2){                               //Make sure that the user is running the code with a configuration file. 
    std::cerr<<"Syntax: "<<argv[0]<<" <Configuration File>"<<std::endl;
    return -1;
  }

  ArrayPack arp;
  std::cerr<<"Argv"<<argv<<std::endl;
  Parameters params(argv[1]);


//Initialise all of the values for an equilibrium or a transient-style run:
//Open data files, set the changing cellsize arrays, do any needed unit conversions. 
//TODO: add flag for transient style along with additional arrays to be opened  

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);  //Text file to save outputs of how much is changing and min and max wtd at various times



  textfile<<"Initialise transient"<<std::endl;



 //load in the data files: ksat, mask, e-folding depth, precipitation, topography, and evaporation files. 


//TODO: Load in 'end' values for each parameter, for transient runs
  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "ksat.nc", "value");   //Units of ksat are m/s. 

  params.ncells_x = arp.ksat.width();  //width and height in number of cells in the array
  params.ncells_y = arp.ksat.height();

  arp.land_mask     = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_mask.nc", "value"); //A binary mask that is 1 where there is land and 0 in the ocean
  arp.fdepth        = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_fdepth_calibrated_1000.nc", "value");  //fslope = 100/(1+150*slope), f>2.5 m. Note this is specific to a 30 arcsecond grid! Other grid resolutions should use different constants. 
  arp.precip        = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_precip.nc",   "value");  //Units: m/yr. 
  arp.temp          = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_temp.nc",   "value");  //Units: degress Celsius
  arp.topo          = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_topo.nc",   "value");  //Units: metres
  arp.starting_evap = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_evap.nc",   "value");  //Units: m/yr
  arp.relhum        = LoadData<float>(params.surfdatadir + params.region + params.time_start + "_relhum.nc",   "value");  //Units: proportion from 0 to 1.

  //Initialise arrays that start out empty:

  arp.wtd    = rd::Array2D<float>(arp.topo,0.0f);  //TODO: load in an actual wtd array. This should come from previous model runs and not start out as 0.
  arp.wtd_old = arp.wtd;
  arp.wtd_mid = arp.wtd;
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

    //cells_per_degree = 120, there are this many 30 arc-second pieces in one degree. (or however many pieces per degree the user has)
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


  //  delete[] arp.cellsize_e_w_metres_S;
//    delete[] arp.cellsize_e_w_metres_N;
  //  delete[] arp.latitude_radians;
//TODO: see which, if any, arrays can be cleared after use to free up memory. How to do this? Above didn't seem to work.


  }

  textfile<<"computed distances, areas, and latitudes"<<std::endl;

  //Change undefined cells to 0
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo(i)<=UNDEF){
      arp.topo(i) = 0;
    }

    //Converting to appropriate time step
    arp.precip(i)        *= (params.deltat/(60*60*24*365)) / 10.0;                  //convert to appropriate units for the time step. 
    arp.starting_evap(i) *= (params.deltat/(60*60*24*365));                  //convert to appropriate units for the time step. 
    
  } 


  arp.check();

//get the starting recharge using precip and evap inputs:

  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.evap(i) = arp.starting_evap(i);               //Use a different array for evaporation here because we want to be able to switch back to the original evaporation in cases where there is surface water for a while, but it goes away.
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);         //Recharge is always positive in locations where there is no surface water. 
  }


//Wtd is 0 in the ocean:
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.land_mask(i) == 0){ 
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
    if(arp.land_mask(i) == 0.0f){ 
      label(i) = dh::OCEAN;
      final_label(i) = dh::OCEAN;
    }
  }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(arp.topo, label, final_label, flowdirs);



 int cycles_done = 0;

 float max_total = 0.0;
 float min_total = 0.0;
 float total_wtd_change = 0.0;
 arp.threshold_array.resize (1000);               //cell area (metres squared)
 for(int i=0;i<1000;i++)
   arp.threshold_array[i] = 100.0;
 int location_in_threshold_array = 0;

  float first_half_thresh = 0.0;
  float second_half_thresh = 0.0;
  float threshold_change = 100.0;
  float wtd_mid_change = 0.0;
  float GW_wtd_change = 0.0;

while(true){

  textfile<<"Cycles done: "<<cycles_done<<std::endl;
  if(cycles_done == 500000)  //How many times to repeat the process. TODO: get a better way of checking for equilibrium. 
    break;

  if((cycles_done % 10) == 0){
    textfile<<"saving partway result"<<std::endl;  
    SaveAsNetCDF(arp.wtd,params.outfilename,"value");

  }


  arp.wtd_old = arp.wtd;
  arp.wtd_mid = arp.wtd;
 //Run the groundwater code to move water

//TODO: Update data arrays every x number of years so they change from starting to end values through time, if we are doing a transient run
  transient(params,arp);

arp.wtd_mid = arp.wtd;

//Move surface water.
  dh::FillSpillMerge(arp.topo, label, final_label, flowdirs, deps, arp.wtd,arp);

//check to see where there is surface water, and adjust how evaporation works at these locations. 
  evaporation_update(params,arp);


  total_wtd_change = 0.0;
  wtd_mid_change = 0.0;
  GW_wtd_change = 0.0;
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    total_wtd_change += fabs(arp.wtd(x,y) - arp.wtd_old(x,y));
    wtd_mid_change += fabs(arp.wtd(x,y) - arp.wtd_mid(x,y));
    GW_wtd_change += fabs(arp.wtd_mid(x,y) - arp.wtd_old(x,y));
   // std::cout<<"total wtd change was "<<total_wtd_change<<" with wtd "<<arp.wtd(x,y)<<" and wtd_old "<<arp.wtd_old(x,y)<<std::endl;
  }

  arp.threshold_array[location_in_threshold_array] = total_wtd_change;
  location_in_threshold_array++;

  if(location_in_threshold_array % 1000 == 0)
    location_in_threshold_array = 0;


  first_half_thresh = 0.0;
  second_half_thresh = 0.0;
  for(int i=0;i<500;i++)
     first_half_thresh += arp.threshold_array[i];
  for(int i=500;i<1000;i++)
    second_half_thresh += arp.threshold_array[i];

  threshold_change = second_half_thresh - first_half_thresh;

    textfile<<"cycles_done "<<cycles_done<<" threshold_change "<<threshold_change<<std::endl;
    textfile<<"first half thresh "<<first_half_thresh<<" second half thresh "<<second_half_thresh<<std::endl;
    textfile<<"total wtd change was "<<total_wtd_change<<" change in GW only was "<<GW_wtd_change<<" and change in SW only was "<<wtd_mid_change<<std::endl;




  arp.wtd_old = arp.wtd;


  cycles_done += 1;
}

//We are finished, save the result.   
  textfile<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,params.outfilename,"value");

  textfile.close();


  return 0;

}
