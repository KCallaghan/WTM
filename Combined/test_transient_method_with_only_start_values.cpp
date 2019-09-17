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

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;


//Initialise all of the values for an transient-style run:
//Open data files, set the changing cellsize arrays, do any needed unit conversions. 
void InitialiseTransient(Parameters &params, ArrayPack &arp){
  std::cout<<"Initialise transient"<<std::endl;

  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");
  arp.land = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); //TODO: for some reason land is loading in with 0s everywhere. Data has both 0s and 1s. Check data export?

  params.width  = arp.ksat.width();
  params.height = arp.ksat.height();


  //Determine area of cells at each latitude
  arp.xlat.resize      (params.height);   // the latitude of each row of cells
  arp.alpha.resize     (params.height);

  const double dy    = 6370000.*6370000.*M_PI/(180.*params.dltxy); //radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.


  //Changing area of cell depending on its latitude. 
  for(unsigned int j=0;j<arp.xlat.size();j++){
    arp.xlat[j]       = (float(j)/params.dltxy+params.sedge)*M_PI/180.;
    //dltxy = 120, there are this many 30 arc-second pieces in one degree. TODO: Change this so that the user can choose a value if they have a different cell size. 
    // j/dltxy gives the number of degrees up from the southern edge, add sedge since the southern edge may not be at 0 latitude. *pi/180 to convert to radians. 
    //xlat is now the latitude in radians. 

    const double xs   = (float(2*j-1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //Latitude number 1 - the southern edge of the cell. TODO: What about when j is 0? fine as long as we're not at the poles.
    const double xn   = (float(2*j+1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //latitude number 2 - the northern edge of the cell. TODO: What about when j is xlat.size()?
    const double area = dy*(std::sin(xn)-std::sin(xs));              //final cell area for that latitude

    // dy    = 6370000.*6370000.*M_PI/(180.*dltxy); //radius of the earth (metres) squared * pi / number of total cells around the world in the y-direction. The part of the area calculation that is constant
  
    //TODO: make sure user-input deltat is also for the correct time step e.g. monthly/annual/daily
   
    arp.alpha[j]      = (params.deltat)/area;  //deltat is the number of seconds per timestep.  
  }

  arp.fdepth   = LoadData<float>(params.surfdatadir + params.time_start   + "_fslope_rotated.nc", "value");
  arp.precip     = LoadData<float>(params.surfdatadir + params.time_start   + "_rech_rotated.nc",   "value");
  arp.temp     = LoadData<float>(params.surfdatadir + params.time_start   + "_temp_rotated.nc",   "value");
  arp.topo     = LoadData<float>(params.surfdatadir + params.time_start   + "_topo_rotated.nc",   "value");


  arp.wtd    = rd::Array2D<float>(arp.topo,0.0f);

  arp.rech    = rd::Array2D<float>(arp.topo,0);
  arp.qtotal    = rd::Array2D<float>(arp.topo,0);

  arp.evap   = LoadData<float>(params.surfdatadir + params.region + "_pretend_evap.nc",   "value");

   arp.temp_east_west.resize(params.height);
  arp.temp_north.resize(params.height);
  arp.temp_south.resize(params.height);




  //Change undefined cells to 0
  for(unsigned int i=0;i<arp.topo.size();i++){
    if(arp.topo(i)<=UNDEF){
      arp.topo(i) = 0;
    }
    //Converting to monthly
    arp.rech(i) /= 12;           //TODO: This should be converted to whatever the timestep is, not necessarily monthly. 
    
    //do we need to do a unit conversion (m -> mm) for rech? For equilibrium we set values >10000 to 0, do we need to do this?

  } 


  arp.check();

    #pragma omp parallel for
  for(unsigned int i=0;i<arp.topo.size();i++){
    arp.evap(i) = arp.evap(i)/100.0f;
    arp.rech(i) = arp.precip(i) - arp.evap(i);
    arp.rech(i) = std::max(arp.rech(i),0.0f);


  }



    for(int y=1;y<params.height-1;y++)     {                      //TODO: how to process edge cells? We can't check all neighbours. Should we check just available neighbours?
      arp.temp_east_west[y] = std::cos(arp.xlat[y]);  
      arp.temp_north[y] = std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.));
      arp.temp_south[y] = std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.));
}


}



//Mini-function that gives the water table height (same as land surface if a wtd isn't loaded in), which = head

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

  for(int y=398;y<params.height-1;y++)
  for(int x=415;x<params.width-1; x++){
    if(arp.ksat(x,y) == 0)
      continue;


    const auto my_head = arp.topo(x,y) + arp.wtd(x,y) + arp.rech(x,y);                     //By doing the head and the kcell like this, we are modifying the wtd one cell at a time and using the modified
    const auto headN   = arp.topo(x,y+1) + arp.wtd(x,y+1) + arp.rech(x,y+1);                       //array for the remainder of the array calculations. We can't do it this way for equilibrium calculations the current 
    const auto headS   = arp.topo(x,y-1) + arp.wtd(x,y-1) + arp.rech(x,y-1);                      //way it's set up, but it makes more sense to me for transient. However, it's easy to switch to calculating and then
    const auto headW   = arp.topo(x-1,y) + arp.wtd(x-1,y) + arp.rech(x-1,y);                       //updating the whole array at once. TODO: Confirm which method is better to use here. 
    const auto headE   = arp.topo(x+1,y) + arp.wtd(x+1,y) + arp.rech(x+1,y);               //Andy says do whole array at once. 

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);


    //TODO: Check all of your signs! I'm not totally sure if it should be - or + here!
    //It seems like we are getting the total which will be discharged from each cell
    double qN = 0;
    double qS = 0;

    double qE = 0;
    double qW = 0;

    arp.qtotal(x,y) = 0;

    //Use equation S3 from the Fan paper. 
    //We get the q in each direction and then add them all together to get the total. This tells us a total of how much water wants to flow into/out of this cell. 



//TODO: consider creating staggered grid with head gradients & kcell averages, and then doing cell-by-cell q discharges. May be faster as we don't calculate each gradient twice. 

//TODO: Separate cell size into a separate calculation. Array of different cell sizes? Similar to the area array



    qN += ((kcellN+my_kcell)/2.)*(headN-my_head) * arp.temp_north[y] * arp.alpha[y]; //North
    qS += ((kcellS+my_kcell)/2.)*(headS-my_head) * arp.temp_south[y] * arp.alpha[y]; //South
    //The extra part is to get to the northern and southern edges of the cell for the width of cell you're touching. 
    qE += (((kcellW+my_kcell)/2.)*(headW-my_head) / arp.temp_east_west[y]) * arp.alpha[y];                             //West
    qW += (((kcellE+my_kcell)/2.)*(headE-my_head) / arp.temp_east_west[y]) * arp.alpha[y];                             //East
   //No extra part needed as the distance doesn't change in this direction, and the length between cells is measured from the middle of the cell. 


//TODO: call this dh or something other than q once we multiply by alpha. 
 
    arp.qtotal(x,y) = (qN + qS + qE + qW)/10.;

 //   std::cout<<"x "<<x<<" y "<<y<<" qtotal "<<arp.qtotal(x,y)<<std::endl;
 //   std::cout<<"x "<<x<<" y "<<y<<" my_kcell "<<my_kcell<<" kcellN "<<kcellN<<" kcellS "<<kcellS<<" kcellE "<<kcellE<<" kcellW "<<kcellW<<std::endl;
//std::cout<<"x "<<x<<" y "<<y<<" my_head "<<my_head<<" headN "<<headN<<" headS "<<headS<<" headE "<<headE<<" headW "<<headW<<std::endl;
//std::cout<<"x "<<x<<" y "<<y<<" topo "<<arp.topo(x,y)<<" wtd "<<arp.wtd(x,y)<<" rech "<<arp.rech(x,y)<<std::endl;

  }



 for(int y=1;y<params.height-1;y++)
 for(int x=1;x<params.width-1; x++){
    if(arp.ksat(x,y) == 0)
      continue;

    if(arp.wtd(x,y)> max_total)
      max_total = arp.wtd(x,y);
    else if(arp.wtd(x,y)< min_total)
      min_total =arp.wtd(x,y);

    arp.wtd(x,y) = arp.wtd(x,y) + arp.qtotal(x,y);   //TODO: For now I am doing it this way but this is something very important to check. When trying to do it with wtdnew, for some reason the output file doesn't want to save, so will need to figure that out if it needs to be done that way. 
 

   total_changes += abs(arp.qtotal(x,y));


  }

std::cout<<"total changes were "<<total_changes<<std::endl;
std::cout<<"max wtd was "<<max_total<<" and min wtd was "<<min_total<<std::endl;

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

    if(iter>=1000)//params.maxiter)// || abs(total_changes) < 20000.0)
      break;

    std::cerr<<"Iteration #: "<<iter<<std::endl;

    cells_left = TransientRun(params, arp, iter,total_changes);
 
    iter++;
  }

  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled-transient.nc","value");


  return 0;
}
