#include "irf.cpp"

#include <richdem/common/timer.hpp>

void initialise(Parameters &params, ArrayPack &arp){
  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);
  //Text file to save outputs of how much is changing and
  //min and max wtd at various times

  if(params.run_type=="transient"){
    textfile<<"Initialise transient"<<std::endl;
    InitialiseTransient(params,arp);
    //compute changing cell size and distances between cells as
    //these change with latitude:
    cell_size_area(params,arp);
    textfile<<"computed distances, areas, and latitudes"<<std::endl;

    //finalise some setup for runoff, labels, etc that is the
    //same for both run types.
    InitialiseBoth(params,arp);
  }
  else if(params.run_type == "equilibrium"){
    textfile<<"Initialise equilibrium"<<std::endl;
    InitialiseEquilibrium(params,arp);
    //compute changing cell size and distances between cells as
    //these change with latitude:
    cell_size_area(params,arp);
    textfile<<"computed distances, areas, and latitudes"<<std::endl;

    //finalise some setup for runoff, labels, etc that is the
    //same for both run types.
    InitialiseBoth(params,arp);
  }
  else if(params.run_type == "test"){
    textfile<<"Initialise test"<<std::endl;
    InitialiseTest(params,arp);
    //compute changing cell size and distances between cells as
    //these change with latitude:
    cell_size_area(params,arp);
    textfile<<"computed distances, areas, and latitudes"<<std::endl;
  }
  else{
    throw std::runtime_error("That was not a recognised run type! \
      Please choose transient or equilibrium.");
  }

  arp.check();
  textfile.close();
}



template<class elev_t>
void update(
  Parameters &params,
  ArrayPack &arp,
  richdem::dephier::DepressionHierarchy<elev_t> &deps
){
  richdem::Timer timer_overall;
  timer_overall.start();

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);

  if(params.run_type == "transient"){
    UpdateTransientArrays(params,arp);
    //linear interpolation of input data from start to end times.
    auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>\
    (arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);
    //with transient runs, we have to redo the depression hierarchy every time,
    //since the topography is changing.
  }

  textfile<<"Cycles done: "<<params.cycles_done<<std::endl;

  //TODO: How should equilibrium know when to exit?
 // if((params.cycles_done % 3) == 0){
    textfile<<"saving partway result."<<std::endl;
    //arp.wtd.printAll();
    string cycles_str = to_string(params.cycles_done);
    arp.wtd.saveGDAL(params.outfilename + cycles_str +".tif");

    //Save the output every 100 iterations, under a new filename
    //so we can compare how the water table has changed through time.
 // }

  arp.wtd_old = arp.wtd;  //These are used to see how much change occurs
  arp.wtd_mid = arp.wtd;  //in FSM vs in the groundwater portion.


  //////////////////////
  // Move groundwater //
  //////////////////////

  time_t now = time(0);
  char* dt = ctime(&now);

  std::cerr << "Before GW time: " << dt << std::endl;

  richdem::Timer time_groundwater;
  time_groundwater.start();


  // AW: Let's rename these to be inner time steps instead of just iterations
  // AW: unless this IS what Kerry intended...
  int iter_count = 0;
  while(iter_count++ < params.maxiter){
    // Apply recharge to the water-table depth grid (wtd)
    // Its clone (wtd_T) is used and updated in the Picard iteration
    #pragma omp parallel for collapse(2)
    for(int y=1;y<params.ncells_y-1;y++)
    for(int x=1;x<params.ncells_x-1; x++){
      if(arp.land_mask(x,y) == 0)          //skip ocean cells
        continue;
      arp.wtd(x,y) = add_recharge(params.deltat, arp.rech(x,y), arp.wtd(x,y), arp.land_mask(x,y), arp.porosity(x,y));
      arp.wtd_T(x,y) = arp.wtd(x,y);
    }
    FanDarcyGroundwater::update(params, arp);
  }


  now = time(0);
  dt = ctime(&now);
  std::cerr << "t GW time = " << time_groundwater.lap() << std::endl;
  std::cerr << "After GW time: " << dt << std::endl;

  arp.wtd_mid = arp.wtd;


  ////////////////////////
  // Move surface water //
  ////////////////////////

  if(params.fsm_on){
    richdem::Timer fsm_timer;
    fsm_timer.start();
    dh::FillSpillMerge(params,deps,arp);

    now = time(0);
    dt = ctime(&now);

    std::cerr << "t FSM time = " << fsm_timer.lap() << std::endl;
    std::cerr << "After FSM time: " << dt << std::endl;
  }


  /////////////////////////////
  // Evaporate surface water //
  /////////////////////////////

  // Check to see where there is surface water, and adjust how evaporation works
  // at these locations.
  richdem::Timer evaporation_timer;
  evaporation_timer.start();
  // evaporation_update(params,arp);
  
  // Evap mode 1: Use the computed open-water evaporation rate
  if(params.evap_mode){
    std::cout<<"updating the evaporation field"<<std::endl;
    #pragma omp parallel for
    for(unsigned int i=0;i<arp.topo.size();i++){
      if(arp.wtd(i)>0)  //if there is surface water present
        arp.rech(i) = arp.precip(i) - arp.open_water_evap(i);
      else{              //water table is below the surface
        arp.rech(i) = arp.precip(i) - arp.starting_evap(i);
        if(arp.rech(i) <0)    //Recharge is always positive.
          arp.rech(i) = 0.0f;
      }
    }
  }
  
  // Evap mode 0: remove all surface water (like Fan Reinfelder et al., 2013)
  else{
    std::cout<<"removing all surface water"<<std::endl;
    #pragma omp parallel for
    for(unsigned int i=0;i<arp.topo.size();i++){
      if(arp.wtd(i)>0)  //if there is surface water present
        arp.wtd(i) = 0;   //use this option when testing GW component alone
      else{              //water table is below the surface
        arp.rech(i) = arp.precip(i) - arp.starting_evap(i);
        if(arp.rech(i) <0)    //Recharge is always positive.
          arp.rech(i) = 0.0f;
      }
    }
  }

  std::cerr << "t Evaporation time = " << evaporation_timer.lap() << std::endl;
  std::cerr << "After evaporation_update: " << dt << std::endl;

  //Print values about the change in water table depth to the text file.
  PrintValues(params,arp);

  arp.wtd_old = arp.wtd;
  params.cycles_done += 1;
  std::cerr << "Done time: " << dt << std::endl;
  std::cerr << "t TWSM update time = " << timer_overall.lap() << std::endl;
}



void run(Parameters &params, ArrayPack &arp){
  //Set the initial depression hierarchy.
  //For equilibrium runs, this is the only time this needs to be done.
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>\
  (arp.topo, arp.cell_area, arp.label, arp.final_label, arp.flowdirs);

  while(true){
    update(params,arp,deps);
    //For transient - user set param that I am setting for now
    //at 50 to get 500 years total.
    if(params.cycles_done == params.total_cycles)
      break;
  }
}


//TODO: how to free up memory? My understanding was that this happens
//automatically when the program terminates, but could be something
//we need to do.
void finalise(Parameters &params, ArrayPack &arp){

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);

  textfile<<"done with processing"<<std::endl;
  string cycles_str = to_string(params.cycles_done);
  arp.wtd.saveGDAL(params.outfilename + cycles_str +".tif");
  //save the final answer for water table depth.

  textfile.close();
}



int main(int argc, char **argv){

  if(argc!=2){
  //Make sure that the user is running the code with a configuration file.
    std::cerr<<"Syntax: "<<argv[0]<<" <Configuration File>"<<std::endl;
    return -1;
  }

  ArrayPack arp;
  std::cerr<<"Argv"<<argv<<std::endl;
  Parameters params(argv[1]);

  initialise(params,arp);
  run(params,arp);
  finalise(params, arp);

  return 0;
}
