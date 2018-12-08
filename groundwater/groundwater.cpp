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


int EquilibriumRun(const Parameters &params, ArrayPack &arp, const int iter){
  int cells_left = 0;

  //We're not allowed to change Parameters, so here we make some copies.
  auto deltat = params.deltat;
  auto alpha  = arp.alpha;

  //I have switched to 30000 iterations before switching to monthly
  //processing, because with 50000 it always seemed to be doing nothing for a
  //long time.

  const int ADJ_COARSE = 1;
  const int ADJ_MEDIUM = 2;
  const int ADJ_FINE   = 3;
  const auto adjustment = ADJ_COARSE;

  //Different increments to use to move water table depth:
  double d0 = 0.005;  //Use /12 for monthly iterations
  double d1 = 0.02;  //No need to use /12 any more, we now switch to monthly within the loop.
  double d2 = 0.1;
  double d3 = 0.25;

  double thres;
  if(adjustment==ADJ_COARSE){
    thres = 0.01;
    d0    = 0.005;
  } else if(adjustment==ADJ_MEDIUM){
    thres = 0.005;
    d0    = d0/5;
  } else if(adjustment==ADJ_FINE){
    thres = 0.001;
    d0    = d0/10;
  } else {
    throw std::runtime_error("Unrecognised increment value.");
  }

  //Do these only once, since we will affect the arrays in a lasting way
  if(iter==30000){
    for(auto i=arp.rech.i0();i<arp.rech.size();i++)
      arp.rech(i) /=12;
    arp.done_old.setAll(false); //We changed the threshold, so we want to recheck all the cells.
  }

  //Do this at each iteration so we don't want to modify the parameter pack
  if(iter>=30000){  //Here we automatically switch to monthly processing,
    std::cerr<<"30000 iterations: adjusting the values for monthly processing."<<std::endl;
    thres  = thres/12.;
    d0     = d0/12.;
    d1     = d1/12.;
    d2     = d2/12.;
    d3     = d3/12.;
    deltat = params.deltat/12.;
    alpha  = arp.alphamonth;
  }

  //!######################################################################################################## change for lakes
  for(auto i=arp.wtd.i0();i<arp.wtd.size();i++)
    arp.wtd(i) = std::min(arp.wtd(i),(float)0); //Any water table above ground gets reset to 0 - this needs to be changed for lakes

  //empty spots in a matrix are zero. (TODO: This note was in the code, not sure what it refers to)


  //TODO: head for ocean cells is empty ie zero, which is correct. I wonder if it is correct for kcell also to be 0?



  for(int y=0;y<params.height;y++)
  for(int x=0;x<params.width; x++){
    if(!arp.land(x,y) || !arp.done_old(x,y))
      continue;

    const auto my_head = head(x,  y,   arp);
    const auto headN   = head(x,  y+1, arp);
    const auto headS   = head(x,  y-1, arp);
    const auto headW   = head(x-1,y,   arp);
    const auto headE   = head(x+1,y,   arp);

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);

    //It seems like we are getting the total which will be discharged from each cell
    double q = 0;
    q += (kcellN+my_kcell)*(headN-my_head) * std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.)); //North   !soo... we're just adding to the total q each time? we're getting a total discharge but not actually moving it in these directions?
    q += (kcellS+my_kcell)*(headS-my_head) * std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.)); //South
    q += (kcellW+my_kcell)*(headW-my_head) / std::cos(arp.xlat[y]);                             //West
    q += (kcellE+my_kcell)*(headE-my_head) / std::cos(arp.xlat[y]);                             //East

    //And we multiply it with alpha, which somehow brings in the timestep?
    q *= arp.alpha[y];
    //I think multiplying it with alpha gets the total that will be discharged as it builds up over that whole time.

    //mm -> m
    // total=rechmean(x,y)*1.e-3 + qlat(x,y)
    //the layer is in metres already.
    const double total = arp.rech(x,y) + q;

   //       As recharge is fixed, the following applies:
   //       (a) if total <0, meaning too much lateral flows i.e., water table is too high.
   //       (b) if total >0, meaning too little lateral flow, i.e., water table is too low.

    //adjustment size depending on how far off it is
    //we use d0-d3 as different size increments of adjustment
    if     (total<-1                                 ) arp.wtd(x,y) += -d3;
    else if(total<-0.25                              ) arp.wtd(x,y) += -d2;
    else if(total<-0.05                              ) arp.wtd(x,y) += -d1;
    else if(total<-thres                             ) arp.wtd(x,y) += -d0;
    //##############################################################################################################change for lakes
    //for now, we have to prevent more from pooling on the surface or we never get equilibrium. BUT this clearly has to change with lakes!
    else if(total>1.    && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d3;
    else if(total>0.25  && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d2;
    else if(total>0.05  && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d1;
    else if(total>thres && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d0;

    if( !(-thres<=total && total<=thres)){  //Did we adjust this cell?
      const int done_n      = arp.done_new(x+1,y)+arp.done_new(x-1,y)+arp.done_new(x,y+1)+arp.done_new(x,y-1)+arp.done_new(x,y);
      cells_left           += done_n;  //If these cells thought they were equilibrated, they thought wrong
      arp.done_new(x+1,y  ) = false;
      arp.done_new(x-1,y  ) = false;
      arp.done_new(x,  y+1) = false;
      arp.done_new(x,  y-1) = false;
      arp.done_new(x,  y  ) = false;
    }
  }

  arp.done_old = arp.done_new;

  return cells_left;
}



int TransientRun(const Parameters &params, ArrayPack &arp, const int iter){
  if(params.interpolated && iter%120==0){
    for(auto i=arp.topo_start.i0();i<arp.topo_start.size();i++){
      arp.topo  (i) = (arp.topo_start  (i) * (1- iter/params.maxiter)) + (arp.topo_end  (i) * (iter/params.maxiter));
      arp.rech  (i) = (arp.rech_start  (i) * (1- iter/params.maxiter)) + (arp.rech_end  (i) * (iter/params.maxiter));
      arp.fdepth(i) = (arp.fslope_start(i) * (1- iter/params.maxiter)) + (arp.fslope_end(i) * (iter/params.maxiter));
      arp.temp  (i) = (arp.temp_start  (i) * (1- iter/params.maxiter)) + (arp.temp_end  (i) * (iter/params.maxiter));
    }

    for(auto i=arp.fdepth.i0();i<arp.fdepth.size();i++){
      if (arp.temp(i)>-5)
        arp.fdepth(i) = arp.fdepth(i);
      else if(arp.temp(i)<-14)
        arp.fdepth(i) = arp.fdepth(i) * (0.17+0.005*arp.temp(i));
      else
        arp.fdepth(i) = arp.fdepth(i) * (1.5 + 0.1*arp.temp(i));
    }
  }


  //!   if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes

  for(int y=0;y<params.height;y++)
  for(int x=0;x<params.width; x++){
    if(!arp.land(x,y))
      continue;

    const auto my_head = head(x,  y,   arp);
    const auto headN   = head(x,  y+1, arp);
    const auto headS   = head(x,  y-1, arp);
    const auto headW   = head(x-1,y,   arp);
    const auto headE   = head(x+1,y,   arp);

    const auto my_kcell = kcell(x,  y,   arp);
    const auto kcellN   = kcell(x,  y+1, arp);
    const auto kcellS   = kcell(x,  y-1, arp);
    const auto kcellW   = kcell(x-1,y,   arp);
    const auto kcellE   = kcell(x+1,y,   arp);

    //It seems like we are getting the total which will be discharged from each cell
    double qnorth = (kcellN+my_kcell)*(headN-my_head) * std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.)); //North
    double qsouth = (kcellS+my_kcell)*(headS-my_head) * std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.)); //South
    double qwest  = (kcellW+my_kcell)*(headW-my_head) / std::cos(arp.xlat[y]);                             //West
    double qeast  = (kcellE+my_kcell)*(headE-my_head) / std::cos(arp.xlat[y]);                             //East

    qnorth *= arp.alpha[y];
    qsouth *= arp.alpha[y];
    qeast  *= arp.alpha[y];
    qwest  *= arp.alpha[y];

    arp.wtd_new(x,y) = arp.wtd(x,y) + qnorth+qsouth+qeast+qwest;    //TODO: Check all of your signs! I'm not totally sure if it should be - or + here!
  }

  arp.wtd = arp.wtd_new;

  return 1;
}



void InitializeCommonBefore(Parameters &params, ArrayPack &arp){
  arp.ksat = LoadData<float>(params.surfdatadir                 + "_ksat.nc", "value");
  arp.land = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value");

  params.width  = arp.ksat.width();
  params.height = arp.ksat.height();

  arp.done_new.resize(arp.topo_start.width(),arp.topo_start.height(),false); //Indicates which cells must still be processed
  arp.done_old.resize(arp.topo_start.width(),arp.topo_start.height(),false); //Indicates which cells must still be processed

  //Determine area of cells at each latitude
  arp.xlat.resize      (arp.topo_start.height());
  arp.alpha.resize     (arp.topo_start.height());
  arp.alphamonth.resize(arp.topo_start.height());

  //Changing area of cell depending on its latitude. TODO: Not totally sure what each of the different variables here represents...
  for(unsigned int j=0;j<arp.xlat.size();j++){
    arp.xlat[j]       = (float(j-2)/params.dltxy+params.sedge)*M_PI/180.;
    const double xs   = (float(2*(j-2)-1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //Latitude number 1
    const double xn   = (float(2*(j-2)+1)/(params.dltxy*2.)+params.sedge)*M_PI/180.; //latitude number 2
    const double area = params.dy*6370000.*(std::sin(xn)-std::sin(xs));              //final cell area for that latitude

    arp.alpha[j]      = 0.5*params.deltat/area;
    arp.alphamonth[j] = 0.5*(params.deltat/12)/area;
  }
}

void InitializeEquilibrium(Parameters &params, ArrayPack &arp){
  arp.fdepth = LoadData<float>(params.surfdatadir + "_fslope_rotated.nc", "value");
  arp.rech   = LoadData<float>(params.surfdatadir + "_rech_rotated.nc",   "value");
  arp.temp   = LoadData<float>(params.surfdatadir + "_temp_rotated.nc",   "value");
  arp.topo   = LoadData<float>(params.surfdatadir + "_topo_rotated.nc",   "value");
  arp.wtd    = rd::Array2D<float>(arp.ksat,0);

  for(auto i=arp.fdepth.i0();i<arp.fdepth.size();i++){
    if(arp.temp(i)>-5)
      arp.fdepth(i) = arp.fdepth(i);
    else if(arp.temp(i)<-14)
      arp.fdepth(i) = arp.fdepth(i) * (0.17+0.005*arp.temp(i));
    else
      arp.fdepth(i) = arp.fdepth(i) * (1.5 + 0.1*arp.temp(i));

    if(arp.fdepth(i)<0.0001)
      arp.fdepth(i) = 0.0001; //Change undefined cells to 0. TODO: Shouldn't there be a cleaner way to do this?
  } 

  for(auto i=arp.rech.i0();arp.rech.size();i++){
    arp.rech(i) *= 1e-3;
    if(arp.rech(i)>=10000)
      arp.rech(i) = 0;
  }  

  //Change undefined cells to 0
  for(auto i=arp.topo.i0();i<arp.topo.size();i++)
    if(arp.topo(i)<=UNDEF)
      arp.topo(i) = 0;

  for(auto i=arp.rech.i0();arp.rech.size();i++) 
    arp.rech(i) = std::max(arp.rech(i), (float)0.);
}

void InitializeTransient(Parameters &params, ArrayPack &arp){
  arp.fslope_start = LoadData<float>(params.surfdatadir + params.time_start + "_fslope_rotated.nc", "value");
  arp.rech_start   = LoadData<float>(params.surfdatadir + params.time_start + "_rech_rotated.nc",   "value");
  arp.temp_start   = LoadData<float>(params.surfdatadir + params.time_start + "_temp_rotated.nc",   "value");
  arp.topo_start   = LoadData<float>(params.surfdatadir + params.time_start + "_topo_rotated.nc",   "value");

  arp.fslope_end   = LoadData<float>(params.surfdatadir + params.time_end   + "_fslope_rotated.nc", "value");
  arp.rech_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_rech_rotated.nc",   "value");
  arp.temp_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_temp_rotated.nc",   "value");
  arp.topo_end     = LoadData<float>(params.surfdatadir + params.time_end   + "_topo_rotated.nc",   "value");

  arp.wtd = LoadData<float>(params.initdatadir + params.region + "_wtd.nc", "value");

  //Change undefined cells to 0
  for(auto i=arp.topo_start.i0();i<arp.topo_start.size();i++)
    if(arp.topo_start(i)<=UNDEF)
      arp.topo_start(i) = 0;
  for(auto i=arp.topo_end.i0();i<arp.topo_end.size();i++)
    if(arp.topo_end(i)<=UNDEF)
      arp.topo_end(i) = 0;  

  //Setting it so recharge can only be positive
  for(auto i=arp.rech_start.i0();arp.rech_start.size();i++) arp.rech_start(i) = std::max(arp.rech_start(i), (float)0.);
  for(auto i=arp.rech_end.i0();  arp.rech_end.size();  i++) arp.rech_end  (i) = std::max(arp.rech_end  (i), (float)0.);    
  
  //Converting to monthly
  for(auto i=arp.rech_start.i0();arp.rech_start.size();i++){
    arp.rech_start(i) /= 12;           
    arp.rech_end  (i) /= 12;
  }        
}

void InitializeCommonAfter(Parameters &params, ArrayPack &arp){
  arp.check();
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

  InitializeCommonBefore(params, arp);

  if(params.run_type=="equilibrium"){
    InitializeEquilibrium(params, arp);
  } else if(params.run_type=="transient"){
    InitializeTransient(params, arp);
  } else {
    throw std::runtime_error("Unrecognised run type!");
  }

  InitializeCommonAfter(params, arp);



  ///////////////////////////////
  //Execution Section

  int iter                  = 0;                           //Number of iterations made
  int cells_to_equilibriate = params.width*params.height;  //Cells left that need to be equilibriated

  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium.
  while(true){
    int cells_left;

    if(params.run_type=="equilibrium" && cells_to_equilibriate<=0.01*params.width*params.height)
      break;
    if(iter>=params.maxiter)
      break;

    std::cerr<<"Iteration #: "<<iter<<std::endl;

    if(params.run_type=="equilibrium")
      cells_left = EquilibriumRun(params, arp, iter);
    else if(params.run_type=="transient")
      cells_left = TransientRun(params, arp, iter);

    if(iter%10000==0){
      //Save the data: wtd
    }

    //TODO
    //call MPI_ALLREDUCE(numbercount,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr) !adds together results from multiple threads

    iter++;
  }

  return 0;
}
