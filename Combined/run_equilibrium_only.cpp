#include "test_equilibrium_only_test_all_at_once_binned.hpp"
//#include "surface_water.hpp"

#include "fill_spill_merge.hpp"

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

namespace rd = richdem;
//namespace dh = richdem::dephier;




int main(int argc, char **argv){

  ArrayPack arp;
  std::cerr<<"Argv"<<argv<<std::endl;
  Parameters params(argv[1]);


  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");
  arp.land_mask = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); 

  params.ncells_x = arp.ksat.width();
  params.ncells_y = arp.ksat.height();

  //Load in data arrays:
  arp.fdepth = LoadData<float>(params.surfdatadir + params.region + "_fslope_rotated.nc", "value");  //fslope = 100/(1+150*slope), f>2.5 m. Note this is specific to a 30 arcsecond grid! Other grid resolutions should use different constants. 
  //TODO: Note the previsor that f>2.5 m! I don't think I have done this in the past, check! 
  arp.precip   = LoadData<float>(params.surfdatadir + params.region + "_rech_rotated.nc",   "value");
  arp.temp   = LoadData<float>(params.surfdatadir + params.region + "_temp_rotated.nc",   "value");
  arp.topo   = LoadData<float>(params.surfdatadir + params.region + "_topo_rotated.nc",   "value");
  arp.evap   = LoadData<float>(params.surfdatadir + params.region + "_pretend_evap.nc",   "value");


  arp.wtd    = rd::Array2D<float>(arp.topo,0.5);
  arp.rech    = rd::Array2D<float>(arp.topo,0);

  arp.head = rd::Array2D<float>(arp.topo,0);        //Just to initialise these - we'll add the appropriate values later. 
  arp.kcell = rd::Array2D<float>(arp.topo,0);

  arp.wtd_change_total = rd::Array2D<float>(arp.topo,0);


  arp.done_new.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed
  arp.done_old.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed


  std::cout<<"loaded all"<<std::endl;

  //Determine area of cells at each latitude

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

  }

  std::cout<<"computed distances, areas, and latitudes"<<std::endl;


  arp.check();


//TODO: USe actual evap layers to get the appropriate precip, this is a placeholder. 
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.evap(i) = arp.evap(i)/100.0f;
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);

  }


//Wtd is 0 in the ocean:
  #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo.isNoData(i)){//} || arp.topo(i)==dh::OCEAN_LEVEL){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      arp.wtd  (i) = 0;
    }
  }

 
  int cells_left = params.ncells_x*params.ncells_y;  //Cells left that need to be equilibriated

 //Run the equilibrium code to move water
  arp.wtd = equilibrium(params,arp, cells_left);

//We are at equilibrium, save the result. 
  
  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled-like-original-more-bins.nc","value");


  return 0;

}
