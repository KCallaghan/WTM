#include "netcdf.hpp"
#include "ArrayPack.hpp"
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

const int ADJ_COARSE = 1;
const int ADJ_MEDIUM = 2;
const int ADJ_FINE   = 3;

const double SNAN = std::numeric_limits<double>::signaling_NaN();

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;



//Mini-function that gives the water table height (same as land surface if a wtd isn't loaded in), which = head
double head(const int x, const int y, const f2d &topo, const f2d &wtd){  
  return topo(x,y) + wtd(x,y);
}

double kcell(const int x, const int y, const f2d &fdepth, const f2d &wtd, const f2d &ksat){
  if(fdepth(x,y)>0){
    if(wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      return fdepth(x,y) * ksat(x,y) * std::exp((wtd(x,y)+1.5)/fdepth(x,y)); //Equation S6 from the Fan paper
    else
      return ksat(x,y) * (wtd(x,y)+1.5+fdepth(x,y));                         //Equation S4 from the Fan paper
  } else {
    return 0;
  }
}


void EquilibriumRun(){
  arp.done_old = arp.done_new;

  //I have switched to 30000 iterations before switching to monthly
  //processing, because with 50000 it always seemed to be doing nothing for a
  //long time.

  const auto adjustment = ADJ_COARSE;

  //Different increments to use to move water table depth:
  double d0 = 0.005;  //Use /12 for monthly iterations
  double d1 = 0.02;  //No need to use /12 any more, we now switch to monthly within the loop.
  double d2 = 0.1;
  double d3 = 0.25;

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


  if(iter==30000){  //Here we automatically switch to monthly processing,
    std::cerr<<"30000 iterations: adjusting the values for monthly processing."<<std::endl;
    thres       = thres/12.;
    d0          = d0/12.;
    d1          = d1/12.;
    d2          = d2/12.;
    d3          = d3/12.;
    deltat      = deltat/12.;
    alpha       = alphamonth;
    rechmean    = rech_month;
    maskold     = 1;            //Because we changed the threshold, so we want to recheck all the cells.
    numbertotal = ntotal;
  }

  //!######################################################################################################## change for lakes
  for(auto i=wtd.i0();i<wtd.size();i++)
    wtd(i) = std::min(wtd(i),0); //Any water table above ground gets reset to 0 - this needs to be changed for lakes

  //empty spots in a matrix are zero. (TODO: This note was in the code, not sure what it refers to)


  //TODO: head for ocean cells is empty ie zero, which is correct. I wonder if it is correct for kcell also to be 0?



  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(!landmask(x,y) || !done_old(x,y))
      continue;

    const auto my_head = head(x,  y,   topo, wtd);
    const auto headN   = head(x,  y+1, topo, wtd);
    const auto headS   = head(x,  y-1, topo, wtd);
    const auto headW   = head(x-1,y,   topo, wtd);
    const auto headE   = head(x+1,y,   topo, wtd);

    const auto my_kcell = kcell(x,y);
    const auto kcellN   = kcell(x,  y+1, fdepth, wtd, ksat);
    const auto kcellS   = kcell(x,  y-1, fdepth, wtd, ksat);
    const auto kcellW   = kcell(x-1,y,   fdepth, wtd, ksat);
    const auto kcellE   = kcell(x+1,y,   fdepth, wtd, ksat);

    //It seems like we are getting the total which will be discharged from each cell
    double q = 0;
    q += (kcellN+my_kcell)*(headN-my_head) * std::cos(xlat(y)+M_PI/(180.*dltxy*2.)) //North   !soo... we're just adding to the total q each time? we're getting a total discharge but not actually moving it in these directions?
    q += (kcellS+my_kcell)*(headS-my_head) * std::cos(xlat(y)-M_PI/(180.*dltxy*2.)) //South
    q += (kcellW+my_kcell)*(headW-my_head) / std::cos(xlat(y))                      //West
    q += (kcellE+my_kcell)*(headE-my_head) / std::cos(xlat(y))                      //East

    //And we multiply it with alpha, which somehow brings in the timestep?
    q *= alpha(y)*q;
    //I think multiplying it with alpha gets the total that will be discharged as it builds up over that whole time.

    //mm -> m
    // total=rechmean(x,y)*1.e-3 + qlat(x,y)
    //the layer is in metres already.
    const double total = rechmean(x,y) + q;

   //       As recharge is fixed, the following applies:
   //       (a) if total <0, meaning too much lateral flows i.e., water table is too high.
   //       (b) if total >0, meaning too little lateral flow, i.e., water table is too low.

    //adjustment size depending on how far off it is
    //we use d0-d3 as different size increments of adjustment
    if     (total<-1                      ) wtd(i,j) += -d3;
    else if(total<-0.25                   ) wtd(i,j) += -d2;
    else if(total<-0.05                   ) wtd(i,j) += -d1;
    else if(total<-thres                  ) wtd(i,j) += -d0;
    //##############################################################################################################change for lakes
    //for now, we have to prevent more from pooling on the surface or we never get equilibrium. BUT this clearly has to change with lakes!
    else if(total>1.    && wtd(i,j)<wtdmax) wtd(i,j) +=  d3;
    else if(total>0.25  && wtd(i,j)<wtdmax) wtd(i,j) +=  d2;
    else if(total>0.05  && wtd(i,j)<wtdmax) wtd(i,j) +=  d1;
    else if(total>thres && wtd(i,j)<wtdmax) wtd(i,j) +=  d0;

    if(-thres<=total && total<=thres){  //Did we adjust this cell?
      total_cells_to_equilibriate--;    //If we didn't adjust the cell, we don't need to think about it again
    } else {
      const int eq_count = arp.done_new(i+1,j)+arp.done_new(i-1,j)+arp.done_new(i,j+1)+arp.done_new(i,j-1)+arp.done_new(i,j);
      total_cells_to_equilibriate -= eq_count;  //If these cells thought they were equilibrated, they thought wrong
      arp.done_new(i+1,j  ) = false;
      arp.done_new(i-1,j  ) = false;
      arp.done_new(i,  j+1) = false;
      arp.done_new(i,  j-1) = false;
      arp.done_new(i,  j  ) = false;
    }
  }

 //maskold = mask !to record which cells still need to be processed.
}



void TransientRun(){
  while(iter++<iterations){
    if(iter%120==0){
      for(auto i=topo_start.i0();i<topo_start.size();i++){
        topo_now  (i) = (topo_start  (i) * (1- iter/iterations)) + (topo_end  (i) * (iter/iterations));
        rech_now  (i) = (rech_start  (i) * (1- iter/iterations)) + (rech_end  (i) * (iter/iterations));
        fslope_now(i) = (fslope_start(i) * (1- iter/iterations)) + (fslope_end(i) * (iter/iterations));
        temp_now  (i) = (temp_start  (i) * (1- iter/iterations)) + (temp_end  (i) * (iter/iterations));
      }

      for(auto i=fdepth.i0();i<fdepth.size();i++){
        if (temp_now(i)>-5)
          fdepth(i) = fslope_now(i);
        else if(temp_now(i)<-14)
          fdepth(i) = fslope_now(i) * (0.17+0.005*temp_now(i));
        else
          fdepth(i) = fslope_now(i) * (1.5 + 0.1*temp_now(i));
      }
    }


    //!   if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes

    const auto my_head = head(x,  y,   topo, wtd);
    const auto headN   = head(x,  y+1, topo, wtd);
    const auto headS   = head(x,  y-1, topo, wtd);
    const auto headW   = head(x-1,y,   topo, wtd);
    const auto headE   = head(x+1,y,   topo, wtd);

    const auto my_kcell = kcell(x,y);
    const auto kcellN   = kcell(x,  y+1, fdepth, wtd, ksat);
    const auto kcellS   = kcell(x,  y-1, fdepth, wtd, ksat);
    const auto kcellW   = kcell(x-1,y,   fdepth, wtd, ksat);
    const auto kcellE   = kcell(x+1,y,   fdepth, wtd, ksat);

    for(int y=0;y<topo.height();y++)
    for(int x=0;x<topo.width(); x++){
      if(landmask(x,y)==0)
        continue;

      //It seems like we are getting the total which will be discharged from each cell
      const double qnorth = (kcellN+my_kcell)*(headN-my_head) * cos(xlat(y)+pi/(180.*delta_xy*2.)); //North
      const double qsouth = (kcellS+my_kcell)*(headS-my_head) * cos(xlat(y)-pi/(180.*delta_xy*2.)); //South
      const double qwest  = (kcellW+my_kcell)*(headW-my_head) / cos(xlat(y));                       //West
      const double qeast  = (kcellE+my_kcell)*(headE-my_head) / cos(xlat(y));                       //East

      qnorth *= alpha(y);
      qsouth *= alpha(y);
      qeast  *= alpha(y);
      qwest  *= alpha(y);

      wtdnew(x,y) = wtd(x,y) + (qnorth(x,y)+qsouth(x,y)+qeast(x,y)+qwest(x,y))    //TODO: Check all of your signs! I'm not totally sure if it should be - or + here!

      wtdnew(x,y+1) = wtd(x,y+1) - qnorth;
      wtdnew(x,y-1) = wtd(x,y-1) - qsouth;
      wtdnew(x-1,y) = wtd(x-1,y) - qwest;
      wtdnew(x+1,y) = wtd(x+1,y) - qeast;
    }

}



int main(int argc, char **argv){
  const double UNDEF  = -1.0e7;
  const double wtdmax = 0;             //Water table depth is 0. TODO: This should probably be changed when I start bringing in the lakes.

  double convergence_threshold = SNAN;    //When to stop calculating - can try with different values and see how it goes. This is the coarsest option.


  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  ArrayPack arp;

  //TODO: region and HAD params

  Parameters params(argv[1]); //TODO

  timer_io.start();
  arp.ksat = LoadData<float>(params.surfdatadir                 + "_ksat.nc");
  arp.land = LoadData<float>(params.initdatadir + params.region + "_mask.nc");

  if(params.run_type=="equilibrium")
    arp.wtd = rd::Array2D<float>(ksat,0);
  else
    arp.wtd = LoadData<float>(params.initdatadir + params.region + "_wtd.nc" );

  arp.topo_start   = LoadData<float>(params.surfdatadir + params.time_start + "_topo_rotated.nc"  );
  arp.rech_start   = LoadData<float>(params.surfdatadir + params.time_start + "_rech_rotated.nc"  );
  arp.fslope_start = LoadData<float>(params.surfdatadir + params.time_start + "_fslope_rotated.nc");
  arp.temp_start   = LoadData<float>(params.surfdatadir + params.time_start + "_temp_rotated.nc"  );

  if(topo_end!=UNINIT_STR){
    arp.fslope_end = LoadData<float>(params.surfdatadir + params.time_end + "_fslope_rotated.nc");
    arp.topo_end   = LoadData<float>(params.surfdatadir + params.time_end + "_topo_rotated.nc"  );
    arp.temp_end   = LoadData<float>(params.surfdatadir + params.time_end + "_temp_rotated.nc"  );
    arp.rech_end   = LoadData<float>(params.surfdatadir + params.time_end + "_rech_rotated.nc"  );
  }

  arp.topo   = arp.topo_start;
  arp.rech   = arp.rech_start;
  arp.fdepth = arp.fslope_start;
  arp.temp   = arp.temp_start;

  arp.equilibrated.resize(topo.width(),topo.height(),false); //Indicates which cells must still be processed

  arp.check();
  timer_io.stop();

  for(auto i=arp.topo_start.i0();i<arp.topo_start.size();i++)       //Change undefined cells to 0
    if(arp.topo_start(i)<=UNDEF)
      arp.topo_start(i) = 0;
  for(auto i=arp.topo_end.i0();i<arp.topo_end.size();i++)           //Change undefined cells to 0
    if(arp.topo_end(i)<=UNDEF)
      arp.topo_end(i) = 0;

  //Determine area of cells at each latitude
  arp.xlat.resize      (arp.topo_start.height());
  arp.alpha.resize     (arp.topo_start.height());
  arp.alphamonth.resize(arp.topo_start.height());

  //Changing area of cell depending on its latitude. TODO: Not totally sure what each of the different variables here represents...
  for(unsigned int j=0;j<xlat.size();j++){
    arp.xlat[j]       = (float(j-2)/dltxy+params.sedge)*pi/180.;
    const double xs   = (float(2*(j-2)-1)/(dltxy*2.)+params.sedge)*pi/180.; //Latitude number 1
    const double xn   = (float(2*(j-2)+1)/(dltxy*2.)+params.sedge)*pi/180.; //latitude number 2
    const double area = dy*6370000.*(sin(xn)-sin(xs));                      //final cell area for that latitude

    arp.alpha[j]      = 0.5*deltat/area
    arp.alphamonth[j] = 0.5*(deltat/12)/area
  }


  //Setting it so recharge can only be positive
  for(auto i=rech_start.i0();rech_start.size();i++) rech_start(i) = std::max(rech_start(i), 0.);
  for(auto i=rech_end.i0();  rech_end.size();  i++) rech_end  (i) = std::max(rech_end  (i), 0.);

  if(params.run_type=="equilibrium"){
    rech_start(i) *= 1e-3;
    if(rech_start(i)>=10000)
      rech_start(i) = 0;
  } else if(params.run_type=="transient"){
    for(auto i=rech_start.i0();rech_start.size();i++){
      rech_start(i) = std::max(rech_start(i),0.);
      rech_start(i) /= 12;           //Converting to monthly
      rech_end(i)   = std::max(rech_end(i),0.);
      rech_end(i)   /= 12;
    }
  } else {
    throw std::runtime_error("Unrecognised run type!");
  }

  for(auto i=fdepth.i0();i<fdepth.size();i++){
    if(temp(i)>-5)
      fdepth(i) = fslope_start(i);
    else if(temp(i,j)<-14)
      fdepth(i) = fslope_start(i) * (0.17+0.005*temp(i));
    else
      fdepth(i) = fslope_start(i) * (1.5 + 0.1*temp(i));

    if(fdepth(i)<0.0001)
      fdepth(i) = 0.0001; //Change undefined cells to 0. TODO: Shouldn't there be a cleaner way to do this?
  }


  int iter = 0;                                  //Count the number of iterations. This is used to know when to switch from annual to monthly cycles.
  int total_cells_to_equilibriate = 0.99*ntotal; //We will use this to get numbertotal less than x % land cells as equilibrium condition

  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium.
  while(true){
    if(total_cells_to_equilibriate--==0)
      break;
    if(iter>=params.maxiter)
      break;

    if(rank==0)
      std::cerr<<"Iteration #: "<<iter<<std::endl;

    if(params.run_type=="equilibrium")
      EquilibriumRun();
    else if(params.run_type=="transient")
      TransientRun();

    if(iter%10000==0){
      //Save the data: wtd
    }

    //TODO
    //call MPI_ALLREDUCE(numbercount,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr) !adds together results from multiple threads

    iter++;
  }

  return 0;
}
