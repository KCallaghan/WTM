int main(){
  rd::Array2D<float> topo_start   = LoadData<float>(surfdatadir + time_start + "_topo_rotated.nc");
  rd::Array2D<float> topo_end     = LoadData<float>(surfdatadir + time_end   + "_topo_rotated.nc");

  rd::Array2D<float> rech_start   = LoadData<float>(surfdatadir + time_start + "_rech_rotated.nc");
  rd::Array2D<float> rech_end     = LoadData<float>(surfdatadir + time_end   + "_rech_rotated.nc");

  rd::Array2D<float> fslope_start = LoadData<float>(surfdatadir + time_start + "_fslope_rotated.nc");
  rd::Array2D<float> fslope_end   = LoadData<float>(surfdatadir + time_end   + "_fslope_rotated.nc");

  rd::Array2D<float> temp_start   = LoadData<float>(surfdatadir + time_start + "_temp_rotated.nc");
  rd::Array2D<float> temp_end     = LoadData<float>(surfdatadir + time_end   + "_temp_rotated.nc");



  rd::Array2D<bool>   equilibrated(topo.width(),topo.height(),false); //Indicates which cells must still be processed
  rd::Array2D<double> head        (topo.width(),topo.height(),   0);  //Indicates which cells must still be processed
  rd::Array2D<bool>   landmask    (topo.width(),topo.height(),   0);  //TODO: Not initialized in Kerry's code

  for(auto i=topo.i0();i<topo.size();i++){       //Change undefined cells to 0
    if(topo_start(i)<=UNDEF)
      topo_start(i)=0;
    if(topo_end(i)<=UNDEF)
      topo_end(i)=0;
  }


  for(auto i=rech_start_read.i0();rech_start_read.size();i++){
    rech_start_read(i) = std::max(rech_start_read(i),0.);
    rech_start_read(i) /= 12;           //Converting to monthly
    rech_end_read(i)   = std::max(rech_end_read(i),0.);
    rech_end_read(i)   /= 12;
  }


}