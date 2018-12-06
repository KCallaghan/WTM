#include "netcdf.hpp"
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


int main(int argc, char **argv){
  const double UNDEF  = -1.0e7;
  const double SEDGE  = -60;           //Southern-most latitude
  const double wtdmax = 0;             //Water table depth is 0. TODO: This should probably be changed when I start bringing in the lakes. 



  const double deltat = 365*24*3600.0; //Seconds in an annual timestep
  
  const double dltxy  = 120;           //There are 120 30 arc-second pieces in one degree
  const double dy     = 6370000*M_PI/(180*dltxy); //radius of the earth * pi / number of possible cells in the y-direction. This should equal the height of each cell in the N-S direction.
  const double dx     = dy;

  double convergence_threshold = SNAN;    //When to stop calculating - can try with different values and see how it goes. This is the coarsest option.

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






  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <OutGraph>"<<std::endl;
    return -1;
  }

  const std::string region      = "000000";    //Folder name for the 4 files. Limit 8 characters.
  const std::string surfdatadir = "surfdata/"; //Folder names for data, limit 60 characters. Topo, fdepth, and rech go here.
  const std::string HAD         = "HADCM3/";
  const std::string initdatadir = "initdata/"; //Mask and wtd go here

  timer_io.start();
  rd::Array2D<float> topo   = LoadData<float>(surfdatadir       + region + "_topo.nc";
  rd::Array2D<float> fdepth = LoadData<float>(surfdatadir       + region + "_fslope.nc";
  rd::Array2D<float> ksat   = LoadData<float>(surfdatadir                + "_ksat.nc";
  rd::Array2D<float> rech   = LoadData<float>(surfdatadir + HAD + region + "_rech.nc";
  rd::Array2D<float> mask   = LoadData<float>(initdatadir       + region + "_mask.nc";
  rd::Array2D<float> wtd    = LoadData<float>(initdatadir       + region + "_wtd.nc";
  rd::Array2D<float> temp   = LoadData<float>(surfdatadir + HAD + region + "_temp.nc";
  timer_io.stop();

  rd::Array2D<bool>   equilibrated(topo.width(),topo.height(),false); //Indicates which cells must still be processed
  rd::Array2D<double> head        (topo.width(),topo.height(),   0);  //Indicates which cells must still be processed
  rd::Array2D<bool>   landmask    (topo.width(),topo.height(),   0);  //TODO: Not initialized in Kerry's code


  for(auto i=topo.i0();i<topo.size();i++)       //Change undefined cells to 0
    if(topo(i)<=UNDEF)
      topo(i)=0;

  wtd.setAll(0);  //No longer using an input water table! 

  for(auto i=rech.i0();rech.size();i++){
    rech(i) = std::max(rech(i),0.);  //Setting it so recharge can only be positive
    rech(i) *= 1e-3;
    if(rech(i)>=10000)
      rech(i) = 0;
  }


  rd::Array2D<float> fdepth_start(fdepth);
  
  for(auto i=fdepth.i0();i<fdepth.size();i++){
    if(temp(i)>-5)
      fdepth_start(i) = fdepth(i);
    else if(temp(i,j)<-14)
      fdepth_start(i) = fdepth(i) * (0.17+0.005*temp(i));
    else
      fdepth_start(i) = fdepth(i) * (1.5 + 0.1*temp(i));

    if(fdepth_start<0.0001)
      fdepth_start = 0.0001; //Change undefined cells to 0. TODO: Shouldn't there be a cleaner way to do this?
  }












 
  //rech_month_read = (rech_read/12.)  






  //Determine area of cells at each latitude
  // std::vector<double> dvec xlat;
  // std::vector<double> dvec alpha;
  // std::vector<double> dvec alphamonth;
  // if(pid .gt. 0) then
  //     nmax = nend(pid) - nini(pid) +1
  //     allocate(xlat(nmax))
  //     allocate(alpha(nmax))
  //     allocate(alphamonth(nmax))

  //     do j=1,nmax  !changing area of cell depending on its latitude. Not totally sure what each of the different variables here represents...
  //         xlat(j) = (float(j+nini(pid)-2)/dltxy+SEDGE)*pi/180.
  //         xs      = (float(2*(j+nini(pid)-2)-1)/(dltxy*2.)+SEDGE)*pi/180. !Latitude number 1
  //         xn      = (float(2*(j+nini(pid)-2)+1)/(dltxy*2.)+SEDGE)*pi/180. !latitude number 2
  //         area    = dy*6370000.*(sin(xn)-sin(xs)) !final cell area for that latitude

  //         alpha(j)      = 0.5*deltat/area 
  //         alphamonth(j) = 0.5*(deltat/12)/area 
  //     end do
  // end if




 
  int iter = 0;            //Count the number of iterations. This is used to know when to switch from annual to monthly cycles. 
  
  int total_cells_to_equilibriate = 0.99*ntotal; //We will use this to get numbertotal less than x % land cells as equilibrium condition
  
  //Start of the main loop with condition - either number of iterations or % equilibrium. I am going for 99% equilibrium. 
  while(total_cells_to_equilibriate-->0 && iter<400000){
    if(rank==0) 
      std::cerr<<"Iteration #: "<<iter<<std::endl;
        // write (15,*) 'PID = 0. Numbertotal:',numbertotal,'ntotal',ntotal

    //I have switched to 30000 iterations before switching to monthly
    //processing, because with 50000 it always seemed to be doing nothing for a
    //long time.

    if(iter==30000){  //Hhere we automatically switch to monthly processing, 
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

    if(iter%10000==0){
      //Save the data: wtd
    }


    numbercount = 0

!######################################################################################################## change for lakes
  for(auto i=wtd.i0();i<wtd.size();i++)
    wtd(i) = std::min(wtd(i),0); //Any water table above ground gets reset to 0 - this needs to be changed for lakes


  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    //equilibrated is used to stop us from repeating all of the steps on cells which have already reached equilibrium. 
    if(fdepth(i,j)>0 && (equilibrated(x,y) || y==1 || y==nmax+1)){ //TODO: What are the y's doing here?
      head(x,y) = topo(x,y) + wtd(x,y); //Gives the water table height (same as land surface if a wtd isn't loaded in), which = head
      if(wtd(x,y)<-1.5) //Work out hydraulic conductivity for each cell
        kcell(x,y) = fdepth(x,y) * ksat(x,y) * std::exp((wtd(x,y)+1.5)/fdepth(x,y)) //Equation S6 from the Fan paper
      else
        kcell(x,y) = ksat(x,y)*(wtd(x,y)+1.5+fdepth(x,y))                           //Equation S4 from the Fan paper 
    }
    //TODO: What about kcell here?
  }


  //empty spots in a matrix are zero. (TODO: This note was in the code, not sure what it refers to)


  //TODO: head for ocean cells is empty ie zero, which is correct. I wonder if it is correct for kcell also to be 0?

  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(!landmask(x,y) || !equilibrated(x,y))
      continue;

      //It seems like we are getting the total which will be discharged from each cell
      double q = 0;
      q  = q + (kcell(i,j+1)+kcell(i,j))*(head(i,j+1)-head(i,j)) * cos(xlat(j)+pi/(180.*dltxy*2.)) //North   !soo... we're just adding to the total q each time? we're getting a total discharge but not actually moving it in these directions?
      q  = q + (kcell(i,j-1)+kcell(i,j))*(head(i,j-1)-head(i,j)) * cos(xlat(j)-pi/(180.*dltxy*2.)) //South
      q  = q + (kcell(i-1,j)+kcell(i,j))*(head(i-1,j)-head(i,j)) / cos(xlat(j))                    //West
      q  = q + (kcell(i+1,j)+kcell(i,j))*(head(i+1,j)-head(i,j)) / cos(xlat(j))                    //East

      //And we multiply it with alpha, which somehow brings in the timestep?
      qlat(x,y) = alpha(j)*q;
      //I think multiplying it with alpha gets the total that will be discharged as it builds up over that whole time. 

      //mm -> m
      // total=rechmean(i,j)*1.e-3 + qlat(i,j)
      //the layer is in metres already. 
      const double total = rechmean(i,j) + qlat(i,j);
      numberold          = numbercount

     //       As recharge is fixed, the following applies:
     //       (a) if total <0, meaning too much lateral flows i.e., water table is too high.
     //       (b) if total >0, meaning too little lateral flow, i.e., water table is too low.

      if(total .lt. -1.) then           //adjustment size depending on how far off it is
          wtd(i,j) = wtd(i,j) -d3       //we use d0-d3 as different size increments of adjustment
          numbercount = numbercount + 1 //and count the cell as not yet being in equilibrium. 
      elseif (total .lt. -0.25) then
          wtd(i,j) = wtd(i,j) -d2
          numbercount = numbercount + 1
      elseif (total .lt. -0.05) then
          wtd(i,j) = wtd(i,j) -d1 
          numbercount = numbercount + 1
      elseif (total .lt. -thres) then
          wtd(i,j) = wtd(i,j) -d0
          numbercount = numbercount + 1

!##############################################################################################################change for lakes
      elseif(total .gt. 1. .and.wtd(i,j).lt.wtdmax) then                   !for now, we have to prevent more from pooling on the surface or we never get equilibrium. BUT this clearly has to change with lakes!
          wtd(i,j) = wtd(i,j) +d3
          numbercount = numbercount + 1
      elseif (total .gt. 0.25 .and.wtd(i,j).lt.wtdmax) then
          wtd(i,j) = wtd(i,j) +d2
          numbercount = numbercount + 1
      elseif (total .gt. 0.05 .and.wtd(i,j).lt.wtdmax) then
          wtd(i,j) = wtd(i,j) +d1 
          numbercount = numbercount + 1
      elseif (total .gt. thres .and.wtd(i,j).lt.wtdmax) then
          wtd(i,j) = wtd(i,j) +d0
          numbercount = numbercount + 1
      !things which are between -thres and thres do not change; they are in equilibrium.
      endif


     if (numberold .ne. numbercount) then
         mask(i+1,j) = 1 !tag cells to show they are still not in equilibrium.
         mask(i-1,j) = 1
         mask(i,j+1) = 1
         mask(i,j-1) = 1
         mask(i,j)   = 1
     endif


        endif
    end do
end do



 maskold = mask !to record which cells still need to be processed.

 numbertotal = numbercount !total number of cells still out of equilibrium this iteration
                     


//Sending and receiving the lines on either side of each section to allow for flow across those lines. 

      call MPI_ALLREDUCE(numbercount,numbertotal,1,MPI_INTEGER,mpi_sum,MPI_COMM_WORLD,ierr) !adds together results from multiple threads

    iter++;
  }


  return 0;
}
