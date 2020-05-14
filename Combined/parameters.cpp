#include "parameters.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Dummy initializer
Parameters::Parameters(){}

// Real initializer
Parameters::Parameters(const std::string config_file){
  std::ifstream fin(config_file);

  if(!fin.good())
    throw std::runtime_error("Failed to read config file!");

  std::string line;
  while(std::getline(fin, line)) {
    if(line.empty())
      continue;

    std::stringstream ss(line);
    std::string key;
    ss>>key;


    //Dummy key to make it easier to alphabetize list below
    if     (key=="")                   {}
    else if(key=="cells_per_degree")   ss>>cells_per_degree;
    else if(key=="deltat")             ss>>deltat;
    else if(key=="infiltration_on")    ss>>infiltration_on;
    else if(key=="maxiter")            ss>>maxiter;
    else if(key=="outfilename")        ss>>outfilename;
    else if(key=="region")             ss>>region;
    else if(key=="run_type")           ss>>run_type;
    else if(key=="southern_edge")      ss>>southern_edge;
    else if(key=="surfdatadir")        ss>>surfdatadir;
    else if(key=="textfilename")       ss>>textfilename;
    else if(key=="time_end")           ss>>time_end;
    else if(key=="time_start")         ss>>time_start;
    else if(key=="total_cycles")       ss>>total_cycles;

    else
      throw std::runtime_error("Unrecognised key!");
  }
}

void Parameters::print() const {
  std::cout<<"c cells_per_degree = "<<cells_per_degree <<std::endl;
  std::cout<<"c deltat           = "<<deltat           <<std::endl;
  std::cout<<"c infiltration_on  = "<<infiltration_on  <<std::endl;
  std::cout<<"c maxiter          = "<<maxiter          <<std::endl;
  std::cout<<"c outfilename      = "<<outfilename      <<std::endl;
  std::cout<<"c region           = "<<region           <<std::endl;
  std::cout<<"c run_type         = "<<run_type         <<std::endl;
  std::cout<<"c southern_edge    = "<<southern_edge    <<std::endl;
  std::cout<<"c surfdatadir      = "<<surfdatadir      <<std::endl;
  std::cout<<"c textfilename     = "<<textfilename     <<std::endl;
  std::cout<<"c time_end         = "<<time_end         <<std::endl;
  std::cout<<"c time_start       = "<<time_start       <<std::endl;
  std::cout<<"c total_cycles     = "<<total_cycles     <<std::endl;
  //TODO: Synchronize with structure
}
