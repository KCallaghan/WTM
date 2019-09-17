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
namespace dh = richdem::dephier;




int main(int argc, char **argv){

  ArrayPack arp;
  std::cerr<<"Argv"<<argv<<std::endl;
  Parameters params(argv[1]);



  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");
  arp.land = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); //TODO: for some reason land is loading in with 0s everywhere. Data has both 0s and 1s. Check data export?

  params.width = arp.ksat.width();
  params.height = arp.ksat.height();

  //Determine area of cells at each latitude
  arp.xlat.resize      (params.height);     // the latitude of each row of cells
  arp.alpha.resize     (params.height);
  arp.alphamonth.resize(params.height);
  arp.temp_east_west.resize(params.height);
  arp.temp_north.resize(params.height);
  arp.temp_south.resize(params.height);


 const double dy    = 6370000.*6370000.*M_PI/(180.*params.dltxy); //radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
 //The part of the area calculation that is constant ^
 //dltxy is the number of cells per degree (e.g. for 30 arc-second cells, dltxy = 120). 


  //Changing area of cell depending on its latitude. 
  for(unsigned int j=0;j<arp.xlat.size();j++){
    arp.xlat[j]       = (float(j)/params.dltxy+params.sedge)*M_PI/180.;
    //dltxy represents the number of cells in one degree. For most of my runs, it is 120 since there are this many 30 arc-second pieces in one degree.
    //j/dltxy gives the number of degrees up from the southern edge, add sedge since the southern edge may not be at 0 latitude. *pi/180 to convert to radians. 
    //xlat is now the latitude in radians. 

    const double xs   = (float(2*j-1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //Latitude number 1 - the southern edge of the cell
    const double xn   = (float(2*j+1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //latitude number 2 - the northern edge of the cell
    const double area = dy*(std::sin(xn)-std::sin(xs));              //final cell area for that latitude



    arp.alpha[j]      = params.deltat/area;    //deltat is the number of seconds per timestep. We assume an annual timestep here, although the user can select a different deltat value, we still /12 to convert to monthly at a set point. 
    arp.alphamonth[j] = (params.deltat/12)/area;   //for when we switch to a monthly version. 
  }


  for(int y=1;y<params.height-1;y++)     {                      //TODO: how to process edge cells? We can't check all neighbours. Should we check just available neighbours?
    arp.temp_east_west[y] = std::cos(arp.xlat[y]);  
    arp.temp_north[y] = std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.));
    arp.temp_south[y] = std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.));
}




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

   arp.total = rd::Array2D<float>(arp.topo,0);


//TODO: Why can't I get these to work in the arraypack?
  rd::Array2D<dh::dh_label_t> label   (arp.topo.width(), arp.topo.height(), dh::NO_DEP ); 



  arp.done_new.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed
  arp.done_old.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed


  std::cout<<"loaded all"<<std::endl;


  arp.check();


  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    arp.evap(i) = arp.evap(i)/100.0f;
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);

  }


//Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(arp.topo.isNoData(i) || arp.topo(i)==dh::OCEAN_LEVEL){ //Ocean Level is assumed to be lower than any other cells (even Death Valley)
      arp.wtd  (i) = 0;
    }
  }

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  int cells_left = params.width*params.height;  //Cells left that need to be equilibriated

 
  arp.wtd = equilibrium(params,arp, cells_left);


  
  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled-like-original-more-bins.nc","value");


  return 0;

}
