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
    arp.rech(i) *= 1e-3;                    //Converting mm to m?  TODO: Check in units to see if this is needed. 
    if(arp.rech(i) >=10000)                 //Impossible values
      arp.rech(i) = 0;
    arp.rech(i) = std::max(arp.rech(i),(float)0.);   //Only positive recharge. TODO: Is it best to keep it this way and handle lake evaporation as a completely separate thing? Or do we want all recharge but we only add it if it's positive?
  }

  std::cout<<"made adjustments to rech and topo"<<std::endl;


  arp.check();
}




int EquilibriumRun(const Parameters &params, ArrayPack &arp, const int loops_passed){

  int cells_left = 0;


  double d0 = 0.005;            //TODO: Consider whether these are the best values to use. Should we have more categories?
  double d1 = 0.02;             //Test larger values for d4/d5  etc
  double d2 = 0.1;
  double d3 = 0.25;
  double thres = 0.01;
  

  //We're not allowed to change Parameters, so here we make some copies.

  auto alpha  = arp.alpha;


 //Do these only once, since we will affect the arrays in a lasting way
  if(loops_passed==500){
    std::cout<<"10 loops"<<std::endl;
    for(auto i=arp.rech.i0();i<arp.rech.size();i++)   //Convert the recharge to monthly. 
      arp.rech(i) /=12;
    arp.done_old.setAll(false); //We changed the threshold, so we want to recheck all the cells.
  }

  //Do this at each iteration so we don't want to modify the parameter pack
  //The selection of 30000 iterations is somewhat arbitrary: We want to do year-long iterations for as 
  //long as this is useful, then switch to month-long to enable us to get close to equlibrium more 
  //quickly. The best value to use likely depends on input data. Original global tests were done with 
  //50000 iterations, and seemed to not do much for a long time, so I switched to 30000. The optimal number probably varies depending on input topography.
  if(loops_passed>=500){  //Here we automatically switch to monthly processing,
    std::cerr<<"10 loops: adjusting the values for monthly processing."<<std::endl;
    thres  = thres/12.;
    d0     = d0/12.;
    d1     = d1/12.;
    d2     = d2/12.;
    d3     = d3/12.;
    alpha  = arp.alphamonth;
  }

  //!######################################################################################################## change for lakes
//  for(auto i=arp.wtd.i0();i<arp.wtd.size();i++)
//    arp.wtd(i) = std::min(arp.wtd(i),(float)0); //Any water table above ground gets reset to 0 - this needs to be changed for lakes


//TODO: Is this a good spot to do evaporation from lakes?
//TODO: Would it be better to have evaporation in the surface water part of code? Should it be in its own file?

// Evaporation(params, arp);







  for(int y=0;y<params.height;y++)                           
  for(int x=0;x<params.width;x++){



    arp.head(x,y) = arp.topo(x,y) + arp.wtd(x,y);        //I removed this from being in its own function since with this method, we want to change everything in the array
                                                               //and THEN update the whole array. So here we calculate whole head array, since below we modify the wtd. 
                                                               //May be other possibilities if we modify both current and neighbour wtds below instead of only current? But would require some thought.



    if(arp.wtd(x,y)<-1.5)            //Work out hydraulic conductivity for each cell
      arp.kcell(x,y) = arp.fdepth(x,y) * arp.ksat(x,y) * std::exp((arp.wtd(x,y)+1.5)/arp.fdepth(x,y)); //Equation S6 from the Fan paper.
    else
      arp.kcell(x,y) = arp.ksat(x,y) * (arp.wtd(x,y)+1.5+arp.fdepth(x,y));  //Equation S4 from the Fan paper. 


    arp.done_new(x,y) = true;  //Set everything to done; those that change will be set back to false again. 

 
}

  for(int y=1;y<params.height-1;y++)                           //TODO: how to process edge cells? We can't check all neighbours. Should we check just available neighbours?
  for(int x=1;x<params.width-1;x++){                           //TODO: will these work until the edges with a working land/ocean mask?
    if(arp.ksat(x,y)==0 || arp.done_old(x,y))    //TODO: not sure how to include the mask. Using ksat instead of land for now due to problems with land layer.
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



   
   
    double q = 0;
    //Use equation S3 from the Fan paper. 
    //We get the q in each direction and then add them all together to get the total. This tells us a total of how much water wants to flow into/out of this cell. 
    q += ((kcellN+my_kcell)/2.)*(headN-my_head) * std::cos(arp.xlat[y]+M_PI/(180.*params.dltxy*2.)); //North   
    q += ((kcellS+my_kcell)/2.)*(headS-my_head) * std::cos(arp.xlat[y]-M_PI/(180.*params.dltxy*2.)); //South
    //The extra part is to get to the northern and southern edges of the cell for the width of cell you're touching. 
    q += ((kcellW+my_kcell)/2.)*(headW-my_head) / std::cos(arp.xlat[y]);                             //West
    q += ((kcellE+my_kcell)/2.)*(headE-my_head) / std::cos(arp.xlat[y]);                             //East
    //No extra part needed as the distance doesn't change in this direction, and the length between cells is measured from the middle of the cell. 

    q *= arp.alpha[y];       // recall arp.alpha[j] = 0.5*deltat/area, where deltat is the number of seconds in a timestep. 
    //So: ksat must originally be in m/s units. * head = m^2/s.
    // We multiply it by the total number of seconds we are processing for to get a total, and divide it by the area of each cell. We now have a unitless value 
    //that represents discharge into/out of this cell, corrected for the amount of time and the area of the cell. 
    //TODO: Why does alpha include the 0.5?



    const double total = arp.rech(x,y) + q;   //TODO: I wonder why we add rech here, as opposed to adding it to head and having it be part of the comparison above?
//TODO: move recharge sum for the whole array at once to above q. 
    if     (total<-1                                 ) arp.wtd(x,y) += -d3;  //Adjust the wtd in the cell according to how much it needs to change to get closer to equilibrium. 
    else if(total<-0.25                              ) arp.wtd(x,y) += -d2;  //TODO: I wonder if it would help to have one still larger adjustment? 
    else if(total<-0.05                              ) arp.wtd(x,y) += -d1;
    else if(total<-thres                             ) arp.wtd(x,y) += -d0;
    else if(total>1.    && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d3;  //Here we are only changing them if the wtd is below the land surface. 
    else if(total>0.25  && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d2;  //TODO: I think it makes more sense without that condition; cases where SW is moving, GW will still move beneath the surface.
    else if(total>0.05  && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d1;  //TODO: Should we have a special class for kcell in such a case? Should kcell max out when wtd = 0?
    else if(total>thres && arp.wtd(x,y)<params.wtdmax) arp.wtd(x,y) +=  d0;


    if( !(-thres<=total && total<=thres)){  //Did we adjust this cell?


      cells_left           += 1;         //If these cells thought they were equilibrated, they thought wrong
      arp.done_new(x+1,y  ) = false;
      arp.done_new(x-1,y  ) = false;
      arp.done_new(x,  y+1) = false;
      arp.done_new(x,  y-1) = false;
      arp.done_new(x,  y  ) = false;
    }

  }


  arp.done_old = arp.done_new;  //Next time, we can skip cells that were within the threshold this time. 

  return cells_left;

}




f2d equilibrium(Parameters &params, ArrayPack &arp, int &cells_left, int &loops_passed){

 InitialiseEquilibrium (params, arp);

  
  int iter = 0;  //TODO: iterations must be counted within combined, else it resets each time. 


while(true){

  cells_left = EquilibriumRun(params,arp,loops_passed);  //TODO: how to get this value to combined, so it knows when to exit there?


  std::cerr<<"Iteration #: "<<iter<<" Cells left: "<<cells_left<<std::endl;

  if(cells_left <= 0.01*params.width*params.height)
    break;

  if(iter >= params.maxiter)
    break;

  iter++;

}

 
  return arp.wtd;

}

