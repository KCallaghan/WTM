#include "parameters.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

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
