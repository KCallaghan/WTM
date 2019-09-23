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


typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;


//Initialise all of the values for an equilibrium-style run:
//Open data files, set the changing cellsize arrays, set appropriate fdepth, do any needed unit conversions. 
void InitialiseEquilibrium (Parameters &params, ArrayPack &arp){

//fdepth = the efolding depth, rate of decay of hydraulic conductivity with depth. It differs depending on slope and temperature. 
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

  std::cout<<"changed to the final fdepths"<<std::endl;


  for(auto i=arp.topo.i0();i<arp.topo.size();i++){
    if(arp.topo(i)<=params.UNDEF)
      arp.topo(i) = 0;
 //   arp.rech(i) *= 1e-3;                    //Converting mm to m?  TODO: Check in units to see if this is needed. 
    if(arp.rech(i) >=10000)                 //Impossible values
      arp.rech(i) = 0;
//    arp.rech(i) = std::max(arp.rech(i),(float)0.);   //Only positive recharge. TODO: Is it best to keep it this way and handle lake evaporation as a completely separate thing? Or do we want all recharge but we only add it if it's positive?
  }

  std::cout<<"made adjustments to rech and topo"<<std::endl;


  arp.check();
}


int EquilibriumRun(const Parameters &params, ArrayPack &arp, const int iter){


  int cells_left = 0;

  double d0 = 0.0005;             //Test larger values for d4/d5  etc
  double d1 = 0.0008;
  double d2 = 0.001;
  double d3 = 0.0015;
  double d4 = 0.002;
  double d5 = 0.003;
  double d6 = 0.005;
  double thres = 0.001;
  


 //Do these only once, since we will affect the arrays in a lasting way
  if(iter==10000){  //TODO: What is a good number of iterations to use here? Is it better to just do monthly from the beginning?
    for(auto i=arp.rech.i0();i<arp.rech.size();i++)   //Convert the recharge to monthly. TODO: Should we move this to combined, or should we do it each time here, since we are now resetting rech each time in combined?
      arp.rech(i) /=12;
    arp.done_old.setAll(false); //We changed the threshold, so we want to recheck all the cells.
  }

  //Do this at each iteration so we don't want to modify the parameter pack. TODO: Or is it better to just do the smaller values from the beginning?
  if(iter>=10000){  //Here we automatically switch to monthly processing,
    thres  = thres/12.;
    d0     = d0/12.;
    d1     = d1/12.;
    d2     = d2/12.;
    d3     = d3/12.;
    d4     = d4/12.;
    d5     = d5/12.;
    d6     = d6/12.;

  
  }

  for(int y=0;y<params.ncells_y;y++)                           
  for(int x=0;x<params.ncells_x;x++){


    arp.head(x,y) = arp.topo(x,y) + arp.wtd(x,y) + arp.rech(x,y);      //trying what happens if I add rech here instead. TODO: It makes sense to me to add rech up here, but seems to go slower to equilibrium. Is it okay to add later instead?

    
    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      arp.kcell(x,y) = arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper.
    else if(arp.wtd(x,y)>0)
       arp.kcell(x,y) = arp.ksat(x,y) * (0.0+1.5+arp.fdepth(x,y)); 
    else
      arp.kcell(x,y) = arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));  //Equation S4 from the Fan paper. 


    arp.done_new(x,y) = true;  //Set everything to done; those that change will be set back to false again. 

 

}
    float maxtotal = 0.0;
    float mintotal = 0.0;
 

  for(int y=1;y<params.ncells_y-1;y++)                           //TODO: how to process edge cells? We can't check all neighbours. Should we check just available neighbours?
  for(int x=1;x<params.ncells_x-1;x++){                           //TODO: will these work until the edges with a working land/ocean mask?
    if(arp.land_mask(x,y)==0 || arp.done_old(x,y))    
      continue;                                   //Skip the cell if it's ocean, or if it's in equilibrium from the last iteration. 


    const auto my_head = arp.head(x,y);              //Get the head and hydraulic conductivity values for the target cell and its neighbours
    const auto headN   = arp.head(x,  y+1);
    const auto headS   = arp.head(x,  y-1);
    const auto headW   = arp.head(x-1,y);
    const auto headE   = arp.head(x+1,y);

    const auto my_kcell = arp.kcell(x,  y);
    const auto kcellN   = arp.kcell(x,  y+1);
    const auto kcellS   = arp.kcell(x,  y-1);
    const auto kcellW   = arp.kcell(x-1,y);
    const auto kcellE   = arp.kcell(x+1,y);
   
    double QN = 0;
    double QS = 0;
    double QE = 0;
    double QW = 0;

    double wtd_change_N = 0;
    double wtd_change_S = 0;
    double wtd_change_E = 0;
    double wtd_change_W = 0;
    //Use equation S3 from the Fan paper. 
    //We get the q in each direction and then add them all together to get the total. This tells us a total of how much water wants to flow into/out of this cell. 
    QN = ((kcellN+my_kcell)/2.)*(headN-my_head)*params.deltat/params.cellsize_n_s_metres;
    QS = ((kcellS+my_kcell)/2.)*(headS-my_head)*params.deltat/params.cellsize_n_s_metres ;
    QE = ((kcellW+my_kcell)/2.)*(headW-my_head)*params.deltat/arp.cellsize_e_w_metres[y] ;
    QW = ((kcellE+my_kcell)/2.)*(headE-my_head)*params.deltat/arp.cellsize_e_w_metres[y];
    
    wtd_change_N  = QN / arp.cell_area[y+1];
    wtd_change_S  = QS/ arp.cell_area[y-1];
    wtd_change_E  = QE/ arp.cell_area[y];
    wtd_change_W  = QW / arp.cell_area[y];


    arp.wtd_change_total(x,y) = (wtd_change_N + wtd_change_S + wtd_change_E + wtd_change_W);

//    arp.total(x,y) = arp.rech(x,y) + q;   //TODO: I wonder why we add rech here, as opposed to adding it to head and having it be part of the comparison above?
//TODO: does it work better to add recharge down here? It seems less logical, but looks like it equilibrates faster?


  }

int d6_total = 0;
int d5_total = 0;
int d4_total = 0;
int d3_total = 0;
int d2_total = 0;
int d1_total = 0;
int d0_total = 0;

  for(int y=1;y<params.ncells_y-1;y++)                           //TODO: how to process edge cells? We can't check all neighbours. Should we check just available neighbours?
  for(int x=1;x<params.ncells_x-1;x++){                           //TODO: will these work until the edges with a working land/ocean mask?
    if(arp.land_mask(x,y)==0 || arp.done_old(x,y))    
      continue;                                   //Skip the cell if it's ocean, or if it's in equilibrium from the last iteration. 

 

//TODO: Check these bins to see if we get something that works better. Although, it probably varies depending on topography - what is the best way to handle this?

  if(arp.wtd_change_total(x,y)<-0.05    ){ 
      arp.wtd(x,y) += -d6;
d6_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-0.02    ){ 
      arp.wtd(x,y) += -d5;
d5_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-0.01   ){
     arp.wtd(x,y) += -d4;
d4_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-0.007  ){ 
      arp.wtd(x,y) += -d3;
d3_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-0.004  ){
      arp.wtd(x,y) += -d2;
d2_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-0.002  ){
     arp.wtd(x,y) += -d1;
d1_total +=1;
   }
    else if(arp.wtd_change_total(x,y)<-thres  ){
      arp.wtd(x,y) += -d0;
d0_total +=1;
   }
  
    
    else if(arp.wtd_change_total(x,y)>0.05     ){ 
      arp.wtd(x,y) +=  d6;
d6_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>0.02     ){ 
      arp.wtd(x,y) +=  d5;
d5_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>0.01    ){ 
      arp.wtd(x,y) +=  d4;
d4_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>0.007    ){ 
      arp.wtd(x,y) +=  d3;
d3_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>0.004    ){ 
      arp.wtd(x,y) +=  d2;
d2_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>0.002   ){ 
      arp.wtd(x,y) +=  d1;
d1_total +=1;
   }
    else if(arp.wtd_change_total(x,y)>thres   ){ 
      arp.wtd(x,y) +=  d0;
d0_total +=1;
 }


    if(arp.wtd_change_total(x,y) >=thres || arp.wtd_change_total(x,y) <= -thres){
      cells_left           += 1;         //If these cells thought they were equilibrated, they thought wrong
      arp.done_new(x,  y  ) = false;
      arp.done_new(x+1,y  ) = false;
      arp.done_new(x-1,y  ) = false;
      arp.done_new(x,  y+1) = false;
      arp.done_new(x,  y-1) = false;

    }
   

    if(arp.wtd_change_total(x,y) > maxtotal)
      maxtotal = arp.wtd_change_total(x,y);

    else if(arp.wtd_change_total(x,y) < mintotal)
      mintotal = arp.wtd_change_total(x,y);

}

  arp.done_old = arp.done_new;  //Next time, we can skip cells that were within the threshold this time. 

  //TODO: Move print statements to output to a text file. 
  std::cout<<"the highest total value was "<<maxtotal<<" and the lowest was "<<mintotal<<std::endl;
  std::cout<<"d6 "<<d6_total<<" d5 "<<d5_total<<" d4 "<<d4_total<<" d3 "<<d3_total<<" d2 "<<d2_total<<" d1 "<<d1_total<<" d0 "<<d0_total<<std::endl;

  return cells_left;

}


f2d equilibrium(Parameters &params, ArrayPack &arp, int &cells_left){

 InitialiseEquilibrium (params, arp);

  
  int iter = 0;  //TODO: iterations must be counted within combined, else it resets each time. 


while(true){

  cells_left = EquilibriumRun(params,arp,iter);  //TODO: how to get this value to combined, so it knows when to exit there?

  std::cerr<<"Iteration #: "<<iter<<" Cells left: "<<cells_left<<std::endl;

  if(cells_left <= 0.05*params.ncells_x*params.ncells_y){
    std::cout<<"Achieved 99.5 percent equilibrium"<<std::endl;
      SaveAsNetCDF(arp.wtd,"test-filled-like-original-more-bins-new.nc","value");

    break;
  }

  if(iter >= params.maxiter){
    std::cout<<"reached max iters"<<std::endl;
      SaveAsNetCDF(arp.wtd,"test-filled-like-original-more-bins-new.nc","value");

    break;
  }

     if((iter % 10000) == 0){
      std::cerr<<"Saving a part-way output"<<std::endl;
      SaveAsNetCDF(arp.wtd,"test-filled-like-original-more-bins-new.nc","value");

    }


  iter++;

}
 
  return arp.wtd;

}