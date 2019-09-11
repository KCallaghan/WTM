#ifndef _debugging_utilities_hpp_
#define _debugging_utilities_hpp_

template<class T>
void PrintDEM(const std::string title, const rd::Array2D<T> &arr, const int width=2){
  return;
  std::cerr<<"\n"<<title<<std::endl;
  std::cerr<<std::setw(2)<<" "<<"    ";
  for(int x=0;x<arr.width();x++)
    std::cerr<<std::setw(width)<<x<<" ";
  std::cerr<<"\n"<<std::endl;
  for(int y=0;y<arr.height();y++){
    std::cerr<<std::setw(2)<<y<<"    ";
    for(int x=0;x<arr.width(); x++){
      if (std::is_same<T, flowdir_t>::value)
        std::cerr<<std::setw(width)<<(int)arr(x,y)<<" ";
      else
        std::cerr<<std::setw(width)<<arr(x,y)<<" ";
    }
    std::cerr<<"     "<<std::setw(2)<<y<<std::endl;
  }
  std::cerr<<"\n"<<std::setw(2)<<" "<<"    ";
  for(int x=0;x<arr.width();x++)
    std::cerr<<std::setw(width)<<x<<" ";
  std::cerr<<"\n"<<std::endl;  
}

template<class elev_t>
void PrintDepressionInfo(const DepressionHierarchy<elev_t> &deps){
  return;
  std::cerr<<"\033[91m######################Depression Info\033[39m"<<std::endl;
  std::cerr<<std::setw(20)<<"Depression"<<std::setw(10)<<"Dep Vol"<<std::setw(10)<<"Water Vol"<<std::endl;
  for(unsigned int d=0;d<deps.size();d++)
    std::cerr<<std::setw(20)<<d<<std::setw(10)<<deps.at(d).dep_vol<<std::setw(10)<<deps.at(d).water_vol<<std::endl;
  std::cerr<<std::endl;
}

template<class elev_t>
void PrintCellsAffectedProfile(const std::vector<int> &cells_affected, const elev_t last_elev, const rd::Array2D<elev_t> &topo){
  std::cerr<<"\nCellsAffectedProfile: ";
  for(const auto x: cells_affected)
    std::cerr<<std::setw(2)<<topo(x)<<" ";
  std::cerr<<std::setw(2)<<last_elev<<std::endl;
}

#endif
