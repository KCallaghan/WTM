#include "parameters.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

void Parameters::print() const {
  std::cout<<"c name        = "<<name        <<std::endl;
  std::cout<<"c iterations  = "<<iterations  <<std::endl;
  std::cout<<"c width       = "<<width       <<std::endl;
  std::cout<<"c height      = "<<height      <<std::endl;
  std::cout<<"c start_name  = "<<start_name  <<std::endl;
  std::cout<<"c end_name    = "<<end_name    <<std::endl;
  std::cout<<"c surfdatadir = "<<surfdatadir <<std::endl;
  std::cout<<"c initdatadir = "<<initdatadir <<std::endl;
  std::cout<<"c delta_xy    = "<<delta_xy    <<std::endl;
  std::cout<<"c dy          = "<<dy          <<std::endl;
  std::cout<<"c deltat      = "<<deltat      <<std::endl;
}

Parameters LoadParameters(const std::string config_file){
  Parameters params;

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

    if(key=="iterations")
      ss>>params.iterations;
    else if(key=="name")
      ss>>params.name;
    else if(key=="width")
      ss>>params.width;
    else if(key=="height")
      ss>>params.height;
    else if(key=="start_name")
      ss>>params.start_name;
    else if(key=="end_name")
      ss>>params.end_name;
    else if(key=="surfdatadir")
      ss>>params.surfdatadir;
    else if(key=="initdatadir")
      ss>>params.initdatadir;
    else if(key=="deltat")
      ss>>deltat;
    else
      throw std::runtime_error("Unrecognised key!");
  } 
}
