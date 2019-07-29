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


		//Parameters for evaporation. TODO: Move to parameter file. 
  double a = 1.26;
  double cp = 1.005; //specific heat of air at 300K. TODO: Should a different temperature be used? Units: kJ/kg.K
  double Lv = 2260;  //latent heat of vapourisation of water. Unit: kJ/kg
  double gamma = cp/Lv; //TODO: Or are we supposed to include pressure and molecular weight ratio as well?
  double p = 1; //TODO: probably we should change the pressure with elevation? Or what?
  double Rv = 461; //specific gas constant of water vapour. Units: J/kg/K
  double S = 100; //Representative of area of water body. TODO: This value is a placeholder. How to do this?
  double u = 1; //TODO: u should be the wind speed. This is a placeholder value; need to load in an array representing this. 
  double f = pow((5 * pow(10.0,6.0) / S), 0.05) * (3.6 + 2.5*u); //Wind speed function
  double RELHUM = 20; //TODO: Load in actual relative humidity arrays; this is a placeholder. 




  arp.ksat = LoadData<float>(params.surfdatadir + params.region + "_ksat.nc", "value");
  arp.land = LoadData<float>(params.initdatadir + params.region + "_mask.nc", "value"); //TODO: for some reason land is loading in with 0s everywhere. Data has both 0s and 1s. Check data export?

  params.width = arp.ksat.width();
  params.height = arp.ksat.height();

  //Determine area of cells at each latitude
  arp.xlat.resize      (params.height);     // the latitude of each row of cells
  arp.alpha.resize     (params.height);
  arp.alphamonth.resize(params.height);

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



  arp.fdepth = LoadData<float>(params.surfdatadir + params.region + "_fslope_rotated.nc", "value");  //fslope = 100/(1+150*slope), f>2.5 m. Note this is specific to a 30 arcsecond grid! Other grid resolutions should use different constants. 
  //TODO: Note the previsor that f>2.5 m! I don't think I have done this in the past, check! 
  arp.rech   = LoadData<float>(params.surfdatadir + params.region + "_rech_rotated.nc",   "value");
  arp.temp   = LoadData<float>(params.surfdatadir + params.region + "_temp_rotated.nc",   "value");
  arp.topo   = LoadData<float>(params.surfdatadir + params.region + "_topo_rotated.nc",   "value");
  arp.wtd    = rd::Array2D<float>(arp.topo,0);
  arp.head = rd::Array2D<float>(arp.topo,0);        //Just to initialise these - we'll add the appropriate values later. 
  arp.kcell = rd::Array2D<float>(arp.topo,0);
  arp.delta = rd::Array2D<float>(arp.topo,0);
  arp.e_a = rd::Array2D<float>(arp.topo,0);
  arp.e = rd::Array2D<float>(arp.topo,0);
  arp.evap = rd::Array2D<float>(arp.topo,0);







  arp.done_new.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed
  arp.done_old.resize(arp.topo.width(),arp.topo.height(),false); //Indicates which cells must still be processed


  std::cout<<"loaded all"<<std::endl;

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

  std::cout<<"changed to final fdepths"<<std::endl;


  for(auto i=arp.topo.i0();i<arp.topo.size();i++){
    if(arp.topo(i)<=params.UNDEF)
      arp.topo(i) = 0;
    arp.rech(i) *= 1e-3;                    //Converting mm to m?  TODO: Check in units to see if this is needed. 
    if(arp.rech(i) >=10000)                 //Impossible values
      arp.rech(i) = 0;
    arp.rech(i) = std::max(arp.rech(i),(float)0.);   //Only positive recharge. TODO: Is it best to keep it this way and handle lake evaporation as a completely separate thing? Or do we want all recharge but we only add it if it's positive?
  }

  std::cout<<"made adjustments to rech and topo"<<std::endl;

//Get the delta values for the domain:
    for(auto i=arp.delta.i0();i<arp.delta.size();i++){  
    	arp.delta(i) = (p*Lv)/(arp.temp(i)*arp.temp(i)*Rv);
    	arp.e_a(i) = 611 * exp((Lv/Rv) * ((1/273.15) - (1/arp.temp(i))));
    	arp.e(i) = (RELHUM/100) * arp.e_a(i);
    	arp.evap(i) = (a/(a-1)) * (gamma/(arp.delta(i) + gamma)) * f * (arp.e_a(i) - arp.e(i));

    }





  arp.check();
}




void Evaporation(const Parameters &params, ArrayPack &arp){

	//TODO: If we choose to ignore the lake area part in the wind function, then we can calculate evaporation fields ahead of time and just load in 
  //an evaporation file, not needing this caluclation in the code and needing one less array loaded in as we don't need to load in the wind speed. 
  //Would it be better to do it like this, or is it better having it in the code so the calculation can all be seen etc? 
  //If we do want to neglect lake size, what should we use for the wind function?
  //If we don't want to, does every cell in the lake see higher evap because of the lake size? That seems weird. 
  //How do we deal with sinuousity of lake shapes if we do want to include it?

//Move evaporation to the surface water. Run x years of groundwater, then move surface water, then allow evap to happen. 
  //Use P-ET for cells with no lakes, P-lake evap for cells with lakes, for the next x years. 
  //The first x years have P-ET everywhere. 

//  double a = 1.26;
//  double cp = 1.005; //specific heat of air at 300K. TODO: Should a different temperature be used? Units: kJ/kg.K
//  double Lv = 2260;  //latent heat of vapourisation of water. Unit: kJ/kg
//  double gamma = cp/Lv; //TODO: Or are we supposed to include pressure and molecular weight ratio as well?
//  double p = 1; //TODO: probably we should change the pressure with elevation? Or what?
//  double Rv = 461; //specific gas constant of water vapour. Units: J/kg/K




  for(auto i=arp.wtd.i0();i<arp.wtd.size();i++){  
  	if(arp.wtd(i) > 0){   //This is a lake cell, evaporation must happen. 
  		std::cout<<"Evaporation of "<<arp.evap(i)<<" is going to happen on this cell of wtd "<<arp.wtd(i)<<std::endl;
  		arp.wtd(i) = std::max(arp.wtd(i) - arp.evap(i),0.0f);                      
  		std::cout<<"After evaportation, wtd is "<<arp.wtd(i)<<std::endl;





  	}

}

}




int EquilibriumRun(const Parameters &params, ArrayPack &arp, const int iter){

  int cells_left = 0;


  double d0 = 0.005;            //TODO: Consider whether these are the best values to use. Should we have more categories?
  double d1 = 0.02;             //Test larger values for d4/d5  etc
  double d2 = 0.1;
  double d3 = 0.25;
  double thres = 0.01;
  

  //We're not allowed to change Parameters, so here we make some copies.

  auto alpha  = arp.alpha;


 //Do these only once, since we will affect the arrays in a lasting way
  if(iter==30000){
    std::cout<<"30000 iters"<<std::endl;
    for(auto i=arp.rech.i0();i<arp.rech.size();i++)   //Convert the recharge to monthly. 
      arp.rech(i) /=12;
    arp.done_old.setAll(false); //We changed the threshold, so we want to recheck all the cells.
  }

  //Do this at each iteration so we don't want to modify the parameter pack
  //The selection of 30000 iterations is somewhat arbitrary: We want to do year-long iterations for as 
  //long as this is useful, then switch to month-long to enable us to get close to equlibrium more 
  //quickly. The best value to use likely depends on input data. Original global tests were done with 
  //50000 iterations, and seemed to not do much for a long time, so I switched to 30000. The optimal number probably varies depending on input topography.
  if(iter>=30000){  //Here we automatically switch to monthly processing,
    std::cerr<<"30000 iterations: adjusting the values for monthly processing."<<std::endl;
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




f2d equilibrium(const Parameters &params, ArrayPack &arp){

  
  int iter = 0;


  int cells_left = params.width*params.height;  //Cells left that need to be equilibriated


while(true){

  cells_left = EquilibriumRun(params,arp,iter);


  std::cerr<<"Iteration #: "<<iter<<" Cells left: "<<cells_left<<std::endl;

  if(cells_left <= 0.01*params.width*params.height)
    break;

  if(iter >= params.maxiter)
    break;

  iter++;

}

  std::cout<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,"test-filled.nc","value");

  return arp.wtd;

}

