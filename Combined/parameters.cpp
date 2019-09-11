#include "parameters.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

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

    if     (key=="")                   {}                 //Dummy key to make it easier to alphabetize list below
    else if(key=="deltat")             ss>>deltat;
    else if(key=="dltxy")              ss>>dltxy;
    else if(key=="end_name")           ss>>end_name;
    else if(key=="HAD")                ss>>HAD;
    else if(key=="initdatadir")        ss>>initdatadir;
    else if(key=="interpolated")       ss>>interpolated;
    else if(key=="iterations")         ss>>iterations;
    else if(key=="maxiter")            ss>>maxiter;
    else if(key=="name")               ss>>name;
    else if(key=="region")             ss>>region;
    else if(key=="run_type")           ss>>run_type;
    else if(key=="sedge")              ss>>sedge;
    else if(key=="start_name")         ss>>start_name;
    else if(key=="surfdatadir")        ss>>surfdatadir;
    else if(key=="time_start")         ss>>time_start;
    else if(key=="time_end")           ss>>time_end;
    else if(key=="wtdmax")             ss>>wtdmax;
    else
      throw std::runtime_error("Unrecognised key!");
  } 


}

void Parameters::print() const {
  std::cout<<"c name        = "<<name        <<std::endl;
  std::cout<<"c iterations  = "<<iterations  <<std::endl;
  std::cout<<"c maxiter     = "<<maxiter     <<std::endl;
  std::cout<<"c start_name  = "<<start_name  <<std::endl;
  std::cout<<"c end_name    = "<<end_name    <<std::endl;
  std::cout<<"c surfdatadir = "<<surfdatadir <<std::endl;
  std::cout<<"c initdatadir = "<<initdatadir <<std::endl;
  std::cout<<"c dltxy       = "<<dltxy       <<std::endl;
//  std::cout<<"c dy          = "<<dy          <<std::endl;
  std::cout<<"c deltat      = "<<deltat      <<std::endl;
  //TODO: Synchronize with structure
}
