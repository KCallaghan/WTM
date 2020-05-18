#include "transient_groundwater.hpp"
#include "fill_spill_merge.hpp"
#include "evaporation.hpp"
#include "add_recharge.hpp"

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

FanDarcyGroundwater gw;

const double UNDEF  = -1.0e7;


///This function initialises those arrays that are needed only for transient 
///model runs. This includes both start and end states for slope, precipitation,
///temperature, topography, ET, and relative humidity. We also have a land vs 
///ocean mask for the end time. It also includes the starting water table depth 
///array, a requirement for transient runs.
///We also calculate the e-folding depth here, using temperature and slope. 
///TODO: the e-folding depth uses some calibration constants that are dependent 
///on cell-size. How to deal with this when a user may
///have differing cell size inputs? Should these be user-set values?
void InitialiseTransient(Parameters &params, ArrayPack &arp){

  arp.vert_ksat = LoadData<float>(params.surfdatadir + params.region + \
  "ksat.nc", "value");   //Units of ksat are m/s. 

//width and height in number of cells in the array
  params.ncells_x = arp.vert_ksat.width();  
  params.ncells_y = arp.vert_ksat.height();


  arp.slope_start         = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_slope.nc",  "value");  //Slope as a value from 0 to 1. 
  arp.precip_start        = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_precip.nc", "value");  //Units: m/yr. 
  arp.temp_start          = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_temp.nc",   "value");  //Units: degress Celsius
  arp.ground_temp_start          = LoadData<float>(params.surfdatadir + \
    params.region + params.time_start + "_ground_temp.nc",   "value");  
    //Units: degress Celsius
  arp.topo_start          = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_topo.nc",   "value");  //Units: metres
  arp.starting_evap_start = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_evap.nc",   "value");  //Units: m/yr
  arp.relhum_start        = LoadData<float>(params.surfdatadir + params.region \
  + params.time_start + "_relhum.nc", "value");  //Units: proportion from 0 to 1
  arp.wind_speed_start  = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_wind_speed.nc", "value");  //Units: m/s


  arp.winter_temp_start  = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_winter_temp.nc", "value");  //Units: degrees celsius



  arp.slope_end         = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_slope.nc",  "value");  //Slope as a value from 0 to 1. 
  arp.land_mask         = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_mask.nc",   "value");  //A binary mask that is 1 where 
  //there is land and 0 in the ocean
  arp.precip_end        = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_precip.nc", "value");  //Units: m/yr. 
  arp.temp_end          = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_temp.nc",   "value");  //Units: degress Celsius
  arp.ground_temp_end          = LoadData<float>(params.surfdatadir + \
    params.region + params.time_end + "_ground_temp.nc",   "value");  
    //Units: degress Celsius
  arp.topo_end          = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_topo.nc",   "value");  //Units: metres
  arp.starting_evap_end = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_evap.nc",   "value");  //Units: m/yr
  arp.relhum_end        = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_relhum.nc", "value");  //Units: proportion from 0 to 1.
  arp.wind_speed_end    = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_wind_speed.nc", "value");  //Units: m/s


  arp.winter_temp_end    = LoadData<float>(params.surfdatadir + params.region + \
  params.time_end + "_winter_temp.nc", "value");  //Units: degrees Celsius
 

  //load in the wtd result from the previous time: 
  arp.wtd    = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_wtd.nc", "value");

  //calculate the fdepth (e-folding depth, representing rate of decay of the 
  //hydraulic conductivity with depth) arrays:
  arp.fdepth_start = rd::Array2D<float>(arp.topo_start,0); 
  arp.fdepth_end   = rd::Array2D<float>(arp.topo_start,0); 
  //TODO: allow user to vary these calibration constants depending on their 
  //input cellsize? Or do some kind of auto variation of them? 
  for(unsigned int i=0;i<arp.topo_start.size();i++){
    if(arp.winter_temp_start(i) > -5)  //then fdepth = f from Ying's equation S7. 
      arp.fdepth_start(i) = std::max(5000/(1+150*arp.slope_start(i)),25.0f);  
    else{ //then fdpth = f*fT, Ying's equations S7 and S8. 
      if(arp.winter_temp_start(i) < -14)
        arp.fdepth_start(i) = (std::max(5000/(1+150*arp.slope_start(i)),25.0f))\
         * (std::max(0.05, 0.17 + 0.005 * arp.winter_temp_start(i)));
      else
        arp.fdepth_start(i) = (std::max(5000/(1+150*arp.slope_start(i)),25.0f))\
         * (std::min(1.0, 1.5 + 0.1 * arp.winter_temp_start(i)));
    }
    if(arp.winter_temp_end(i) > -5)  //then fdepth = f from Ying's equation S7. 
      arp.fdepth_end(i) = std::max(5000/(1+150*arp.slope_end(i)),25.0f);  
    else{ //then fdpth = f*fT, Ying's equations S7 and S8. 
      if(arp.winter_temp_end(i) < -14)
        arp.fdepth_end(i) = (std::max(5000/(1+150*arp.slope_end(i)),25.0f)) * \
      (std::max(0.05, 0.17 + 0.005 * arp.winter_temp_end(i)));
      else
        arp.fdepth_end(i) = (std::max(5000/(1+150*arp.slope_end(i)),25.0f)) * \
      (std::min(1.0, 1.5 + 0.1 * arp.winter_temp_end(i)));
    }
  }

//initialise the arrays to be as at the starting time:
  arp.fdepth        = arp.fdepth_start;
  arp.precip        = arp.precip_start;
  arp.temp          = arp.temp_start;
  arp.winter_temp   = arp.winter_temp_start;
  arp.topo          = arp.topo_start;
  arp.starting_evap = arp.starting_evap_start;
  arp.relhum        = arp.relhum_start;
  arp.slope         = arp.slope_start;
  arp.evap          = arp.starting_evap;
  arp.ground_temp   = arp.ground_temp_start;
  arp.wind_speed    = arp.wind_speed_start;
}


///This function initialises those arrays that are needed only for equilibrium 
///model runs. 
///This includes a single array for each of slope, precipitation, temperature, 
///topography, ET, land vs ocean mask, and relative humidity. 
///It also includes setting the starting water table depth array to 
///zero everywhere.
///We also calculate the e-folding depth here, using temperature and slope. 
///TODO: the e-folding depth uses some calibration constants that are dependent 
///on cell-size. How to deal with this when a user may
///have differing cell size inputs? Should these be user-set values?
void InitialiseEquilibrium(Parameters &params, ArrayPack &arp){

  arp.vert_ksat = LoadData<float>(params.surfdatadir + params.region + \
  "ksat.nc", "value");   //Units of ksat are m/s. 

//width and height in number of cells in the array
  params.ncells_x = arp.vert_ksat.width();  
  params.ncells_y = arp.vert_ksat.height();

  arp.slope         = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_slope.nc",  "value");  //Slope as a value from 0 to 1. 
  arp.land_mask     = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_mask.nc",   "value"); 
  //A binary mask that is 1 where there is land and 0 in the ocean
  arp.precip        = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_precip.nc", "value");  //Units: m/yr. 
  arp.temp          = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_temp.nc",   "value");  //Units: degress Celsius
  arp.ground_temp          = LoadData<float>(params.surfdatadir + \
  params.region + params.time_start + "_ground_temp.nc",   "value");  
  //Units: degress Celsius
  arp.topo          = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_topo.nc",   "value");  //Units: metres
  arp.starting_evap = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_evap.nc",   "value");  //Units: m/yr
  arp.relhum        = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_relhum.nc", "value");  //Units: proportion from 0 to 1.
  arp.wind_speed    = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_wind_speed.nc", "value");  //Units: m/s


  arp.winter_temp    = LoadData<float>(params.surfdatadir + params.region + \
  params.time_start + "_winter_temp.nc", "value");  //Units: degrees Celsius



  arp.wtd           = rd::Array2D<float>(arp.topo,0.0);  
  //we start with a water table at the surface for equilibrium runs. 
  arp.evap          = arp.starting_evap;

  arp.fdepth   = rd::Array2D<float>(arp.topo,0); 
  for(unsigned int i=0;i<arp.topo.size();i++){
    //TODO: allow user to vary these calibration constants depending on their 
    //input cellsize? Or do some kind of auto variation of them? 
    if(arp.winter_temp(i) > -5)  //then fdepth = f from Ying's equation S7. 
      arp.fdepth(i) = std::max(1000/(1+150*arp.slope(i)),25.0f);  
    else{ //then fdpth = f*fT, Ying's equations S7 and S8. 
      if(arp.winter_temp(i) < -14)
        arp.fdepth(i) = (std::max(1000/(1+150*arp.slope(i)),25.0f)) * \
      (std::max(0.05, 0.17 + 0.005 * arp.winter_temp(i)));
      else
        arp.fdepth(i) = (std::max(1000/(1+150*arp.slope(i)),25.0f)) * \
      (std::min(1.0, 1.5 + 0.1 * arp.winter_temp(i)));
    }
  }
}



///Here, we use the number of cells per degree (a user-defined value), 
///the southern-most latitude of the domain (also user-defined), 
///and the radius of the Earth to calculate the latitude of each row of cells, 
///the size of a cell in the N-S and E-W directions, and the area of each cell. 
void cell_size_area(Parameters &params, ArrayPack &arp){
 //compute changing cell size and distances between 
  //cells as these change with latitude:

  const float earth_radius = 6371000.; //metres

//distance between lines of latitude is a constant. 
  params.cellsize_n_s_metres = (earth_radius*(M_PI/180.))\
  /params.cells_per_degree; 

//initialise some arrays
  arp.latitude_radians.resize       (params.ncells_y);  
  //the latitude of each row of cells
  arp.cellsize_e_w_metres.resize    (params.ncells_y);  
  //size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cellsize_e_w_metres_N.resize  (params.ncells_y);  
  //size of a cell in the east-west direction at the 
  //northern edge of the cell (metres)
  arp.cellsize_e_w_metres_S.resize  (params.ncells_y);  
  //size of a cell in the east-west direction at the 
  //southern edge of the cell (metres)
  arp.cell_area.resize              (params.ncells_y);  
  //cell area (metres squared)


  for(unsigned int j=0;j<arp.latitude_radians.size();j++){
    //latitude at the centre of a cell:
    arp.latitude_radians[j] = (float(j)/params.cells_per_degree + \
    params.southern_edge)*(M_PI/180.); 
    //southern edge of the domain in degrees, plus the number of cells up 
    //from this location/the number of cells per degree, converted to radians.

    //cells_per_degree = 120, there are this many 30 arc-second pieces in 
    //one degree. (or however many pieces per degree the user has)
    //j/cells_per_degree gives the number  of degrees up from the southern edge,
    // add southern_edge since the southern edge may not be at 0 latitude. 
    //*pi/180 to convert to radians. 
    //latitude_radians is now the latitude in radians. 

    //latitude at the southern edge of a cell (subtract half a cell):
    double latitude_radians_S  = ((float(j) - 0.5)\
      /params.cells_per_degree+params.southern_edge)*(M_PI/180.); 
    //latitude at the northern edge of a cell (add half a cell):
    double latitude_radians_N  = ((float(j) + 0.5)\
      /params.cells_per_degree+params.southern_edge)*(M_PI/180.); 

    //distance between lines of longitude varies with latitude. 
    //This is the distance at the centre of a cell for a given latitude:
    arp.cellsize_e_w_metres[j] = earth_radius* \
    std::cos(arp.latitude_radians[j])*(M_PI/180.)/params.cells_per_degree;

    //distance at the northern edge of the cell for the given latitude:
    arp.cellsize_e_w_metres_N[j] = earth_radius* \
    std::cos(latitude_radians_N)*(M_PI/180.)/params.cells_per_degree;
    //distance at the southern edge of the cell for the given latitude:
    arp.cellsize_e_w_metres_S[j] = earth_radius* \
    std::cos(latitude_radians_S)*(M_PI/180.)/params.cells_per_degree;

    //cell area computed as a trapezoid, using unchanging north-south distance, 
    //and east-west distances at the northern and southern edges of the cell:
    arp.cell_area[j] = params.cellsize_n_s_metres* \
    (arp.cellsize_e_w_metres_N[j] + arp.cellsize_e_w_metres_S[j])/2;


//TODO: see which, if any, arrays can be cleared after use to free up memory. 
    //How to do this?
  }
}


///This function initialises those arrays that are used for both equilibrium 
///and transient model runs. This includes arrays that start off with zero 
///values, as well as the label, final_label, and flowdirs arrays. 
void InitialiseBoth(const Parameters &params, ArrayPack &arp){

  arp.ksat = LoadData<float>(params.surfdatadir + params.region + \
  "horizontal_ksat.nc", "value");   //Units of ksat are m/s. 
  arp.porosity    = LoadData<float>(params.surfdatadir + params.region + \
  "porosity.nc", "value");  //Units: unitless


  //Set arrays that start off with zero or other values, 
  //that are not imported files. Just to initialise these - 
  //we'll add the appropriate values later. 

  //These two are just informational, to see how much change 
  //happens in FSM vs in groundwater
  arp.wtd_old            = arp.wtd;                               
  arp.wtd_mid            = arp.wtd;

  arp.runoff             = rd::Array2D<float>(arp.ksat,0);
  arp.head               = rd::Array2D<float>(arp.ksat,0);        
  arp.kcell              = rd::Array2D<float>(arp.ksat,0);

  //Several arrays that are used for calculations of evaporation
  arp.evap               = rd::Array2D<float>(arp.ksat,0);        
  arp.e_sat              = rd::Array2D<float>(arp.ksat,0);  
  arp.e_a                = rd::Array2D<float>(arp.ksat,0);  
  arp.surface_evap       = rd::Array2D<float>(arp.ksat,0); 

  //These are used to see how much change occurred in infiltration 
  //and updating lakes portions of the code. Just informational.  
  arp.infiltration_array = rd::Array2D<float>(arp.ksat,0);        
  arp.surface_array      = rd::Array2D<float>(arp.ksat,0);

  arp.rech               = rd::Array2D<float>(arp.ksat,0);

    arp.stability_time_seconds               = rd::Array2D<float>(arp.ksat,0);


  //This array is used to store the values of how much the water table will 
  //change in one iteration, then adding it to wtd gets the new wtd.
  arp.wtd_change_total   = rd::Array2D<float>(arp.ksat,0.0f);     

  //These are populated during the calculation of the depression hierarchy:
  arp.label              = rd::Array2D<dh_label_t>   \
  (params.ncells_x, params.ncells_y, dh::NO_DEP ); 
  //No cells are part of a depression
  arp.final_label        = rd::Array2D<dh_label_t>   \
  (params.ncells_x, params.ncells_y, dh::NO_DEP ); 
  //No cells are part of a depression
  arp.flowdirs           = rd::Array2D<rd::flowdir_t>\
  (params.ncells_x, params.ncells_y, rd::NO_FLOW); 
  //No cells flow anywhere

  //Change undefined cells to 0
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo(i)<=UNDEF){
      arp.topo(i) = 0;
    }
    //Converting units to appropriate time step
 //   arp.precip(i)        *= (params.deltat/(60*60*24*365));                 
   // arp.starting_evap(i) *= (params.deltat/(60*60*24*365));                 

  } 

//get the starting runoff using precip and evap inputs:
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.rech(i) = (arp.precip(i)-arp.starting_evap(i));
    if(arp.rech(i) <0)    //Recharge is always positive. 
      arp.rech(i) = 0.0f;
    if(arp.porosity(i) <=0 )
      arp.porosity(i) = 0.0000001; //not sure why it is sometimes processing cells with 0 porosity?
  }

//Wtd is 0 in the ocean:
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.land_mask(i) == 0){ 
      arp.wtd  (i) = 0;
    }
  }

//Label the ocean cells. This is a precondition for 
  //using `GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.label.size();i++){
    if(arp.land_mask(i) == 0.0f){ 
      arp.label(i) = dh::OCEAN;
      arp.final_label(i) = dh::OCEAN;
    }
  }
}


///In transient runs, we adjust the input arrays via a 
//linear interpolation from the start state to the end state at each iteration. 
///We do so here, and also reset the label and flow direction arrays, 
///since the depression hierarchy needs to be 
///recalculated due to the changed topography. 
void UpdateTransientArrays(const Parameters &params, ArrayPack &arp){
  for(unsigned int i=0;i<arp.topo.size();i++){

    arp.fdepth(i)         = (arp.fdepth_start(i)        * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.fdepth_end(i)       \
       * (params.cycles_done/params.total_cycles));
    arp.precip(i)         = (arp.precip_start(i)        * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.precip_end(i)       \
       * (params.cycles_done/params.total_cycles));
    arp.temp(i)           = (arp.temp_start(i)          * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.temp_end(i)         \
       * (params.cycles_done/params.total_cycles));
    arp.winter_temp(i)    = (arp.winter_temp_start(i)   * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.winter_temp_end(i)  \
       * (params.cycles_done/params.total_cycles));

    arp.topo(i)           = (arp.topo_start(i)          * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.topo_end(i)         \
       * (params.cycles_done/params.total_cycles));
    arp.starting_evap(i)  = (arp.starting_evap_start(i) * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.starting_evap_end(i)\
       * (params.cycles_done/params.total_cycles));
    arp.relhum(i)         = (arp.relhum_start(i)        * \
      (1-(params.cycles_done/params.total_cycles))) + (arp.relhum_end(i)       \
       * (params.cycles_done/params.total_cycles));

    //Converting to appropriate time step
  //  arp.precip(i)        *= (params.deltat/(60*60*24*365));                  
   // arp.starting_evap(i) *= (params.deltat/(60*60*24*365));                  
 
    arp.label(i)        = dh::NO_DEP; //No cells are part of a depression
    arp.final_label(i)  = dh::NO_DEP; //No cells are part of a depression
    arp.flowdirs(i)     = rd::NO_FLOW; //No cells flow anywhere
  } 

  #pragma omp parallel for
  for(unsigned int i=0;i<arp.label.size();i++){
    if(arp.land_mask(i) == 0.0f){ 
      arp.label(i) = dh::OCEAN;
      arp.final_label(i) = dh::OCEAN;
    }
  }
}


///In this function, we use a few of the variables that were created for 
///informational purposes to help us understand how much the water table 
///is changing per iteration, and where in 
///the code that change is occurring. We print these values to a text file.
void PrintValues(Parameters &params, ArrayPack &arp){

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);  

  params.abs_total_wtd_change = 0.0;
  params.abs_wtd_mid_change = 0.0;
  params.abs_GW_wtd_change = 0.0;
  params.total_wtd_change = 0.0;
  params.wtd_mid_change = 0.0;
  params.GW_wtd_change = 0.0;

  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    params.abs_total_wtd_change += fabs(arp.wtd(x,y)     - arp.wtd_old(x,y));
    params.abs_GW_wtd_change    += fabs(arp.wtd(x,y)     - arp.wtd_mid(x,y));
    params.abs_wtd_mid_change   += fabs(arp.wtd_mid(x,y) - arp.wtd_old(x,y));
    params.total_wtd_change     += (arp.wtd(x,y)         - arp.wtd_old(x,y));
    params.GW_wtd_change        += (arp.wtd(x,y)         - arp.wtd_mid(x,y));
    params.wtd_mid_change       += (arp.wtd_mid(x,y)     - arp.wtd_old(x,y));
    params.infiltration_change  += arp.infiltration_array(x,y);
    params.surface_change       += arp.surface_array(x,y);
  }

  textfile<<"params.cycles_done "<<params.cycles_done<<std::endl;
  textfile<<"total wtd change was "<<params.total_wtd_change<<\
  " change in GW only was "<<params.GW_wtd_change<<\
  " and change in SW only was "<<params.wtd_mid_change<<std::endl;
  textfile<<"absolute value total wtd change was "<<params.abs_total_wtd_change\
  <<" change in GW only was "<<params.abs_GW_wtd_change<<\
  " and change in SW only was "<<params.abs_wtd_mid_change<<std::endl;
  textfile<<"the change in infiltration was "<<params.infiltration_change\
  <<" and change to surface water was "<<params.surface_change<<std::endl;
  textfile.close();
}

