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

const double UNDEF  = -1.0e7;
float cellsize_n_s_metres = 0;


typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;


//Initialise all of the values for an transient-style run:
//Open data files, set the changing cellsize arrays, do any needed unit conversions. 
void InitialiseTransient(Parameters &params, ArrayPack &arp){
  std::cout<<"Initialise transient"<<std::endl;

 //load in the data files: ksat and the mask file. 

  //TODO: check units of ksat. I think it is m/s. 
  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");

  params.ncells_x  = arp.ksat.width();  //width and height in number of cells in the array
  params.ncells_y = arp.ksat.height();

//Load in data arrays
  arp.land_mask = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); 
  arp.fdepth    = LoadData<float>(params.surfdatadir + params.time_start   + "_fslope_rotated.nc", "value");
  arp.precip    = LoadData<float>(params.surfdatadir + params.time_start   + "_rech_rotated.nc",   "value");
  arp.temp      = LoadData<float>(params.surfdatadir + params.time_start   + "_temp_rotated.nc",   "value");
  arp.topo      = LoadData<float>(params.surfdatadir + params.time_start   + "_topo_rotated.nc",   "value");
  arp.evap      = LoadData<float>(params.surfdatadir + params.region + "_pretend_evap.nc",   "value");
  arp.wtd       = rd::Array2D<float>(arp.topo,0.0f);  //TODO: load in an actual wtd array


  for(auto i=arp.fdepth.i0();i<arp.fdepth.size();i++){  //Equation S8 in the Fan paper
    if(arp.temp(i)<-14) {
      auto fT = (0.17 + 0.005*arp.temp(i));
      fT = std::max(fT,0.05);                       //The equation specifies fT>=0.05.
      arp.fdepth(i) = arp.fdepth(i) * fT;
    }
    else if (arp.temp(i) <= -5){
      auto fT = (1.5 + 0.1*arp.temp(i));
      fT = std::min(fT,1.);                       //The equation specifies fT<=1.
      arp.fdepth(i) = arp.fdepth(i) * fT;
    }

    if(arp.fdepth(i)<0.0001)      //TODO: I believe this shouldn't be necessary once I make the needed changes to the original in f array.
      arp.fdepth(i) = 0.0001;
  }



  //initialise arrays that start out empty:
  arp.rech    = rd::Array2D<float>(arp.topo,0);
  arp.wtd_change_total    = rd::Array2D<float>(arp.topo,0);
  arp.stability_time_seconds = rd::Array2D<float>(arp.topo,0);


  //compute changing cell size and distances between cells as these change with latitude:
  const float earth_radius = 6371000.; //metres

  cellsize_n_s_metres = (earth_radius*(M_PI/180.))/params.cells_per_degree; //distance between lines of latitude is a constant. 


  arp.latitude_radians.resize       (params.ncells_y);   // the latitude of each row of cells
  arp.cellsize_e_w_metres.resize    (params.ncells_y);   //size of a cell in the east-west direction at the centre of the cell (metres)
  arp.cellsize_e_w_metres_N.resize  (params.ncells_y);   //size of a cell in the east-west direction at the northern edge of the cell (metres)
  arp.cellsize_e_w_metres_S.resize  (params.ncells_y);   //size of a cell in the east-west direction at the southern edge of the cell (metres)
  arp.cell_area.resize              (params.ncells_y);   //cell area (metres squared)


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
    arp.cell_area[j] = cellsize_n_s_metres* (arp.cellsize_e_w_metres_N[j] + arp.cellsize_e_w_metres_S[j])/2;

  }

  std::cout<<"computed distances, areas, and latitudes"<<std::endl;

  //Change undefined cells to 0
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo(i)<=UNDEF){
      arp.topo(i) = 0;
    }
    if(arp.ksat(i)<=0){
      arp.ksat(i) = 1E-10; //TODO:What is an appropriate minimum value?
    }
    //Converting to monthly
    arp.rech(i) *= 100;           //TODO: This should be converted to whatever the timestep is, not necessarily monthly. 
    
    //do we need to do a unit conversion (m -> mm) for rech? For equilibrium we set values >10000 to 0, do we need to do this?

  } 

  arp.check();


//get the starting recharge using precip and evap inputs:
  //TODO: USe actual evap layers to get the appropriate precip, this is a placeholder. 

    #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.evap(i) = arp.evap(i)/100.0f;
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);


  }

//We have finished initialising all of the arrays and data that will be used. 

}


//Mini-function that gives the kcell, which changes through time as the water table depth changes:
double kcell(const int x, const int y, const ArrayPack &arp){
  if(arp.fdepth(x,y)>0){
    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      return arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper
    else if(arp.wtd(x,y) > 0)
      return arp.ksat(x,y) * (0+1.5+arp.fdepth(x,y)); 
    else
      return arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));                             //Equation S4 from the Fan paper
  } else {
    return 0;
  }
}




int TransientRun(const Parameters &params, ArrayPack &arp, const int iter, double total_changes){

  std::cout<<"TransientRun"<<std::endl;
 
  total_changes = 0.0;
  float max_total = 0.0;
  float min_total = 0.0;
  float local_fdepth = 0.0;
  float local_kcell = 0.0;
  float local_wtd = 0.0;
  float max_change = 0.0;

  int thres_001 = 0;
  int thres_01 = 0;
  int thres_1 = 0;

  float stability_min = 100000000000000000000.0;
  float stability_max = 0.0;


//cycle through the entire array, calculating how much the water table changes in each cell per iteration. 
  for(int y=1;y<params.ncells_y-1;y++)
  for(int x=1;x<params.ncells_x-1; x++){
    if(arp.land_mask(x,y) == 0)
      continue;


    const auto my_head = arp.topo(x,y) + arp.wtd(x,y) + arp.rech(x,y);                    
    const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1) + arp.rech(x,y+1);              
    const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1) + arp.rech(x,y-1);              
    const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y) + arp.rech(x-1,y);              
    const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y) + arp.rech(x+1,y);              

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);

if(my_kcell == 0)
  std::cout<<"0 kcell with fdepth "<<arp.fdepth(x,y)<<" ksat "<<arp.ksat(x,y)<<" wtd "<<arp.wtd(x,y)<<std::endl;

    arp.stability_time_seconds(x,y) = 0.5 * (arp.cell_area[y]*arp.cell_area[y]) / my_kcell;  //Von Neumann stability analysis to get stable timestep
    //TODO: implement stable timestep


    if(arp.stability_time_seconds(x,y)<stability_min){
      stability_min = arp.stability_time_seconds(x,y);
      local_fdepth = arp.fdepth(x,y);
      local_kcell = my_kcell;
      local_wtd = arp.wtd(x,y);
    }
    if(arp.stability_time_seconds(x,y)>stability_max)
      stability_max = arp.stability_time_seconds(x,y);



    double QN = 0;
    double QS = 0;
    double QE = 0;
    double QW = 0;

    double wtd_change_N = 0;
    double wtd_change_S = 0;
    double wtd_change_E = 0;
    double wtd_change_W = 0;


    arp.wtd_change_total(x,y) = 0;



//TODO: consider creating staggered grid with head gradients & kcell averages, and then doing cell-by-cell q discharges. May be faster as we don't calculate each gradient twice. 


//discharge per unit area in each direction:
    QN = ((kcellN+my_kcell)/2.)*(headN-my_head)*params.deltat/cellsize_n_s_metres;
    QS = ((kcellS+my_kcell)/2.)*(headS-my_head)*params.deltat/cellsize_n_s_metres ;
    QE = ((kcellW+my_kcell)/2.)*(headW-my_head)*params.deltat/arp.cellsize_e_w_metres[y];
    QW = ((kcellE+my_kcell)/2.)*(headE-my_head)*params.deltat/arp.cellsize_e_w_metres[y] ;

//multiply by number of seconds we are moving water for, divide by distance water is travelling, divide by area of cell into which it is flowing: 
    wtd_change_N  = QN / arp.cell_area[y+1];
    wtd_change_S  = QS / arp.cell_area[y-1];
    wtd_change_E  = QE / arp.cell_area[y];
    wtd_change_W  = QW / arp.cell_area[y];

 
    arp.wtd_change_total(x,y) = (wtd_change_N + wtd_change_S + wtd_change_E + wtd_change_W);

//    std::cout<<"x "<<x<<" y "<<y<<" qtotal "<<arp.qtotal(x,y)<<std::endl;
   // std::cout<<"x "<<x<<" y "<<y<<" temp_north "<<arp.temp_north[y]<<" cellsize e w "<<arp.cellsize_e_w_metres[y]<<std::endl;
//std::cout<<"x "<<x<<" y "<<y<<" my_head "<<my_head<<" headN "<<headN<<" headS "<<headS<<" headE "<<headE<<" headW "<<headW<<std::endl;
//std::cout<<"x "<<x<<" y "<<y<<" topo "<<arp.topo(x,y)<<" wtd "<<arp.wtd(x,y)<<" rech "<<arp.rech(x,y)<<std::endl;

    if(arp.wtd(x,y)> max_total)
      max_total = arp.wtd(x,y);
    else if(arp.wtd(x,y)< min_total)
      min_total =arp.wtd(x,y);

    if(abs(arp.wtd_change_total(x,y)) > max_change)
      max_change =abs(arp.wtd_change_total(x,y));
  


  }


 for(int y=1;y<params.ncells_y-1;y++)
 for(int x=1;x<params.ncells_x-1; x++){
    if(arp.ksat(x,y) == 0)
      continue;

    arp.wtd(x,y) = arp.wtd(x,y) + arp.wtd_change_total(x,y);   
 
    total_changes += arp.wtd_change_total(x,y);

  }

std::cout<<"total changes were "<<total_changes<<std::endl;
std::cout<<"max wtd was "<<max_total<<" and min wtd was "<<min_total<<std::endl;
std::cout<<"stability_min "<<stability_min<<" stability_max "<<stability_max<<std::endl;
std::cout<<"min kcell "<<local_kcell<<" min wtd "<<local_wtd<<" min fdepth "<<local_fdepth<<std::endl;

std::cout<<"max change was "<< max_change<<std::endl;


  return 1;
}



int main(int argc, char **argv){
  richdem::Timer timer_io;

  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <Configuration File>"<<std::endl;
    return -1;
  }

  ArrayPack arp;

  //TODO: region and HAD params

  Parameters params(argv[1]);

  ///////////////////////////////
  //Initialization section

  InitialiseTransient(params, arp);


  ///////////////////////////////
  //Execution Section

  int iter                  = 0;                           //Number of iterations made
  double total_changes      = 0.0;
  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium.
  while(true){
    int cells_left;

    if((iter % 100) == 0){
      std::cerr<<"Saving a part-way output"<<std::endl;
      SaveAsNetCDF(arp.wtd,"test-filled-transient-cells-done-100-year.nc","value");

    }

    if(iter>=500000)//params.maxiter)// || abs(total_changes) < 20000.0)
      break;

    std::cerr<<"Iteration #: "<<iter<<std::endl;

    cells_left = TransientRun(params, arp, iter,total_changes);
 
    iter++;
  }

  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled-transient-cells-done-100-year.nc","value");


  return 0;
}
