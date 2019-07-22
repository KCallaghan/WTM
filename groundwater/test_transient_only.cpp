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
  arp.xlat.resize      (arp.topo_start.height());   // the latitude of each row of cells
  arp.alpha.resize     (arp.topo_start.height());

  const double dy    = 6370000.*6370000.*M_PI/(180.*params.dltxy); //radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.


  //Changing area of cell depending on its latitude. TODO: Not totally sure what each of the different variables here represents...
  for(unsigned int j=0;j<arp.xlat.size();j++){
    arp.xlat[j]       = (float(j)/params.dltxy+params.sedge)*M_PI/180.;
    //dltxy = 120, there are this many 30 arc-second pieces in one degree. TODO: Change this so that the user can choose a value if they have a different cell size. 
    // j/dltxy gives the number of degrees up from the southern edge, add sedge since the southern edge may not be at 0 latitude. *pi/180 to convert to radians. 
    //xlat is now the latitude in radians. 

    const double xs   = (float(2*j-1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //Latitude number 1 - the southern edge of the cell. TODO: What about when j is 0?
    const double xn   = (float(2*j+1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //latitude number 2 - the northern edge of the cell. TODO: What about when j is xlat.size()?
    const double area = dy*(std::sin(xn)-std::sin(xs));              //final cell area for that latitude

    // dy    = 6370000.*6370000.*M_PI/(180.*dltxy); //radius of the earth (metres) squared * pi / number of total cells around the world in the y-direction. The part of the area calculation that is constant
  
    //TODO: make sure user-input deltat is also for the correct time step e.g. monthly/annual/daily
   
    arp.alpha[j]      = 0.5*params.deltat/area;  //deltat is the number of seconds per timestep.  
  }

  arp.fslope_start = LoadData<float>(params.surfdatadir + params.time_start + "_fslope_rotated.nc", "value"); //fslope = 100/(1+150*slope), f>2.5 m. Note this is specific to a 30 arcsecond grid! Other grid resolutions should use different constants. 
    //TODO: Note the previsor that f>2.5 m! I don't think I have done this in the past, check! 

  arp.rech_start   = LoadData<float>(params.surfdatadir + params.time_start + "_rech_rotated.nc",   "value");
  arp.temp_start   = LoadData<float>(params.surfdatadir + params.time_start + "_temp_rotated.nc",   "value");
  arp.topo_start   = LoadData<float>(params.surfdatadir + params.time_start + "_topo_rotated.nc",   "value");

  arp.fslope_end   = LoadData<float>(params.surfdatadir + params.time_end   + "_fslope_rotated.nc", "value");
  arp.rech_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_rech_rotated.nc",   "value");
  arp.temp_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_temp_rotated.nc",   "value");
  arp.topo_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_topo_rotated.nc",   "value");


  arp.fdepth   = LoadData<float>(params.surfdatadir + params.time_end   + "_fslope_rotated.nc", "value");
  arp.rech     = LoadData<float>(params.surfdatadir + params.time_end   + "_rech_rotated.nc",   "value");
  arp.temp     = LoadData<float>(params.surfdatadir + params.time_end   + "_temp_rotated.nc",   "value");
  arp.topo     = LoadData<float>(params.surfdatadir + params.time_end   + "_topo_rotated.nc",   "value");


  arp.wtd    = rd::Array2D<float>(arp.topo_start,0); //TODO: Change to an actual wtd file


  //Change undefined cells to 0
  for(auto i=arp.topo_start.i0();i<arp.topo_start.size();i++){
    if(arp.topo_start(i)<=UNDEF){
      arp.topo_start(i) = 0;
      arp.topo_end(i) = 0;  
    }
    //Setting it so recharge can only be positive
    arp.rech_start(i) = std::max(arp.rech_start(i), (float)0.);
    arp.rech_end  (i) = std::max(arp.rech_end  (i), (float)0.);  
    //Converting to monthly
    arp.rech_start(i) /= 12;           //TODO: This should be converted to whatever the timestep is, not necessarily monthly. 
    arp.rech_end  (i) /= 12;
    //do we need to do a unit conversion (m -> mm) for rech? For equilibrium we set values >10000 to 0, do we need to do this?

  } 


  arp.check();

}



//Mini-function that gives the water table height (same as land surface if a wtd isn't loaded in), which = head
double head(const int x, const int y, const ArrayPack &arp){  
  return arp.topo(x,y) + arp.wtd(x,y);
}

double kcell(const int x, const int y, const ArrayPack &arp){
  if(arp.fdepth(x,y)>0){
    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      return arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper
    else
      return arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));                             //Equation S4 from the Fan paper
  } else {
    return 0;
  }
}




int TransientRun(const Parameters &params, ArrayPack &arp, const int iter){

  std::cout<<"TransientRun"<<std::endl;
  if(iter % 120 == 0){              //increment the inputs between the start and end values
    for(auto i=arp.topo_start.i0();i<arp.topo_start.size();i++){
      arp.topo  (i) = (arp.topo_start  (i) * (1- iter/params.maxiter)) + (arp.topo_end  (i) * (iter/params.maxiter));
      arp.rech  (i) = (arp.rech_start  (i) * (1- iter/params.maxiter)) + (arp.rech_end  (i) * (iter/params.maxiter));
      arp.fdepth(i) = (arp.fslope_start(i) * (1- iter/params.maxiter)) + (arp.fslope_end(i) * (iter/params.maxiter));
      arp.temp  (i) = (arp.temp_start  (i) * (1- iter/params.maxiter)) + (arp.temp_end  (i) * (iter/params.maxiter));
    }


//fdepth = the efolding depth, rate of decay of hydraulic conductivity with depth. It differs depending on slope and temperature. 

    for(auto i=arp.fdepth.i0();i<arp.fdepth.size();i++){  //re-calculate the fdepth when the input values get updated 
      if (arp.temp(i)<-14)  {                    //Equation S8 in the Fan paper
      auto fT = (0.17 + 0.005*arp.temp(i));
      fT = std::max(fT,0.05);                       //The equation specifies fT>=0.05.
      arp.fdepth(i) = arp.fdepth(i) * fT;
      }      
      else if(arp.temp(i)<= -5){
      auto fT = (1.5 + 0.1*arp.temp(i));
      fT = std::min(fT,1.);                       //The equation specifies fT<=1.
      arp.fdepth(i) = arp.fdepth(i) * fT;
      }
      if(arp.fdepth(i)<0.0001)      //I believe this shouldn't be necessary once I make the needed changes to the original in f array.
        arp.fdepth(i) = 0.0001;
    }
  }


  //!   if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes

  for(int y=0;y<params.height;y++)
  for(int x=0;x<params.width; x++){
    if(!arp.land(x,y))
      continue;

    const auto my_head = head(x,  y,   arp);                        //By doing the head and the kcell like this, we are modifying the wtd one cell at a time and using the modified
    const auto headN   = head(x,  y+1, arp);                        //array for the remainder of the array calculations. We can't do it this way for equilibrium calculations the current 
    const auto headS   = head(x,  y-1, arp);                        //way it's set up, but it makes more sense to me for transient. However, it's easy to switch to calculating and then
    const auto headW   = head(x-1,y,   arp);                        //updating the whole array at once. TODO: Confirm which method is better to use here. 
    const auto headE   = head(x+1,y,   arp);

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);


    //TODO: Check all of your signs! I'm not totally sure if it should be - or + here!
    //It seems like we are getting the total which will be discharged from each cell
    double q = 0;
    //Use equation S3 from the Fan paper. 
    //We get the q in each direction and then add them all together to get the total. This tells us a total of how much water wants to flow into/out of this cell. 

    q += (kcellN+my_kcell)*(headN-my_head) * std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.)); //North
    q += (kcellS+my_kcell)*(headS-my_head) * std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.)); //South
    //The extra part is to get to the northern and southern edges of the cell for the width of cell you're touching. 
    q += (kcellW+my_kcell)*(headW-my_head) / std::cos(arp.xlat[y]);                             //West
    q += (kcellE+my_kcell)*(headE-my_head) / std::cos(arp.xlat[y]);                             //East
   //No extra part needed as the distance doesn't change in this direction, and the length between cells is measured from the middle of the cell. 
                         

    q *= arp.alpha[y];    // recall arp.alpha[j] = 0.5*deltat/area, where deltat is the number of seconds in a timestep. 
    //So: ksat must originally be in m/s units. * head = m^2/s.
    // We multiply it by the total number of seconds we are processing for to get a total, and divide it by the area of each cell. We now have a unitless value 
    //that represents discharge into/out of this cell, corrected for the amount of time and the area of the cell. 
    //TODO: Why does alpha include the 0.5?



    //TODO: Where is rech getting added in for transient runs?

  //  arp.wtd_new(x,y) = arp.wtd(x,y) + q;   //Hmm, actually, we have a separate wtd array here, so we are updating the whole array at once. TODO: Is this way best, or should we update the array cell-by-cell?
 
    arp.wtd(x,y) = arp.wtd(x,y) + q;   //TODO: For now I am doing it this way but this is something very important to check. When trying to do it with wtdnew, for some reason the output file doesn't want to save, so will need to figure that out if it needs to be done that way. 


  }

  //arp.wtd = arp.wtd_new;

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

  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium.
  while(true){
    int cells_left;

    if(iter>=params.maxiter)
      break;

    std::cerr<<"Iteration #: "<<iter<<std::endl;

    
    cells_left = TransientRun(params, arp, iter);

    if(iter%10000==0){
      //Save the data: wtd
    }

    //TODO
    //call MPI_ALLREDUCE(numbercount,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr) !adds together results from multiple threads

    iter++;
  }

  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled-transient.nc","value");


  return 0;
}
