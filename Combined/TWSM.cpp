#include "irf.cpp"




void initialise(Parameters &params, ArrayPack &arp){
  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);  //Text file to save outputs of how much is changing and min and max wtd at various times

  if(params.run_type=="transient"){
    textfile<<"Initialise transient"<<std::endl;
    InitialiseTransient(params,arp);
  }
  else if(params.run_type == "equilibrium"){
    textfile<<"Initialise equilibrium"<<std::endl;
    InitialiseEquilibrium(params,arp);
  }
  else{
    throw std::runtime_error("That was not a recognised run type! Please choose transient or equilibrium.");
  }

  //compute changing cell size and distances between cells as these change with latitude:
  cell_size_area(params,arp);
  textfile<<"computed distances, areas, and latitudes"<<std::endl;

  //finalise some setup for runoff, labels, etc that is the same for both run types.
  InitialiseBoth(params,arp);

  arp.check();
  textfile.close();
}





template<class elev_t>
void update(Parameters &params, ArrayPack &arp, richdem::dephier::DepressionHierarchy<elev_t>   &deps){

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);  //Text file to save outputs of how much is changing and min and max wtd at various times

  if(params.run_type == "transient"){
    UpdateTransientArrays(params,arp);  //linear interpolation of input data from start to end times. 
    auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(arp, arp.label, arp.final_label, arp.flowdirs);    //with transient runs, we have to redo the depression hierarchy every time, since the topography is changing. 
  }


  textfile<<"Cycles done: "<<params.cycles_done<<std::endl;

//TODO: How should equilibrium know when to exit?

  if((params.cycles_done % 100) == 0){
    textfile<<"saving partway result"<<std::endl;  
    string cycles_str = to_string(params.cycles_done);
    SaveAsNetCDF(arp.wtd,params.outfilename + cycles_str +".nc","value");  //Save the output every 100 iterations, under a new filename so we can compare how the water table has changed through time. 
  }

  arp.wtd_old = arp.wtd;  //These are used to see how much change occurs in FSM vs in the groundwater portion. 
  arp.wtd_mid = arp.wtd;

  //Move surface water
  dh::FillSpillMerge(params,deps,arp);

  arp.wtd_mid = arp.wtd;

  //Run the groundwater code to move water
  groundwater(params,arp);

  //check to see where there is surface water, and adjust how evaporation works at these locations. TODO is this still appropriate with the new method of cell-by-cell evaporation?
  evaporation_update(params,arp);

  //Print values about the change in water table depth to the text file. 
  PrintValues(params,arp);
  
  arp.wtd_old = arp.wtd;

  params.cycles_done += 1;
}






void run(Parameters &params, ArrayPack &arp){
  //Set the initial depression hierarchy. For equilibrium runs, this is the only time this needs to be done. 
  auto deps = dh::GetDepressionHierarchy<float,rd::Topology::D8>(arp, arp.label, arp.final_label, arp.flowdirs);

  while(true){
    update(params,arp,deps);
    if(params.cycles_done == params.total_cycles)  //For transient - user set param that I am setting for now at 50 to get 500 years total. 
      break;
  }
}





//TODO: how to free up memory? My understanding was that this happens automatically when the program terminates, but could be something we need to do. 
void finalise(Parameters &params, ArrayPack &arp){

  ofstream textfile;
  textfile.open (params.textfilename, std::ios_base::app);  //Text file to save outputs of how much is changing and min and max wtd at various times

  textfile<<"done with processing"<<std::endl;  
  SaveAsNetCDF(arp.wtd,params.outfilename,"value");  //save the final answer for water table depth. 

  textfile.close();
}





int main(int argc, char **argv){

  if(argc!=2){                               //Make sure that the user is running the code with a configuration file. 
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