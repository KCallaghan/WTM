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


  while(iter++<iterations){
    //!   if (pid .gt. 0) wtd = min(wtd,0.) !any water table above ground gets reset to 0 - this needs to be changed for lakes

    if(iter%120==0){
      for(int i=0;i<topo_start_read.size();i++){
        topo_now(i)   = (topo_start_read(i)   * (1- iter/iterations)) + (topo_end_read(i)   * (iter/iterations))
        rech_now(i)   = (rech_start_read(i)   * (1- iter/iterations)) + (rech_end_read(i)   * (iter/iterations))
        fslope_now(i) = (fslope_start_read(i) * (1- iter/iterations)) + (fslope_end_read(i) * (iter/iterations))
        temp_now(i)   = (temp_start_read(i)   * (1- iter/iterations)) + (temp_end_read(i)   * (iter/iterations))
      }

      for(int i=0;i<fdepth_now.size();i++){
        if (temp_now(i)>-5)
            fdepth_now(i) = fslope_now(i);
        else if(temp_now(i)<-14)
            fdepth_now(i) = fslope_now(i) * (0.17+0.005*temp_now(i));
        else
            fdepth_now(i) = fslope_now(i) * (1.5 + 0.1*temp_now(i));
      }
    }


    for(int y=0;y<topo.height();y++)
    for(int x=0;x<topo.width();x++){
      if(fdepth_sent<=0)  //Maskold is used to stop us from repeating all of the steps on cells which have already reached equilibrium. 
        continue;
      //wtd(i,j) = wtd(i,j) + rech_sent(i,j)  !adding in any recharge for that time step. 

      head(x,y) = topo_sent(x,y) + wtd(x,y) + rech_sent(x,y); //Gives the water table height (same as land surface if a wtd isn't loaded in), which = head

      if(wtd(x,y)<-1.5)   //Work out hydraulic conductivity for each cell
        kcell(x,y) = fdepth_sent(x,y) *ksat(x,y)*exp((wtd(x,y)+1.5)/fdepth_sent(x,y)) !This is equation S6 from the paper
      else
        kcell(x,y) = ksat(x,y)*(wtd(x,y)+1.5+fdepth_sent(x,y)) !equation S4 from the paper 

      //TODO
      //*******************************************************************************************************************************
      //                  else                                                      !here I am treating the surface water as another layer of groundwater
      //                     kcell(x,y) = 1*(wtd(x,y) + 1.5+fdepth_sent(x,y))                       !this is kind of just a placefiller equation to make it do something above the surface. I'm sure we can do better. 
    }

    for(int y=0;y<topo.height();y++)
    for(int x=0;x<topo.width(); x++){
      if(landmask(x,y)==0)
        continue;

      //It seems like we are getting the total which will be discharged from each cell
      qnorth = (kcell(x,y+1)+kcell(x,y))*(head(x,y+1)-head(x,y)) * cos(xlat(y)+pi/(180.*delta_xy*2.)) //North
      qsouth = (kcell(x,y-1)+kcell(x,y))*(head(x,y-1)-head(x,y)) * cos(xlat(y)-pi/(180.*delta_xy*2.)) //South
      qwest  = (kcell(x-1,y)+kcell(x,y))*(head(x-1,y)-head(x,y)) / cos(xlat(y))                       //West
      qeast  = (kcell(x+1,y)+kcell(x,y))*(head(x+1,y)-head(x,y)) / cos(xlat(y))                       //East

      const double qlat_north = alpha(y)*qnorth;
      const double qlat_south = alpha(y)*qsouth;
      const double qlat_east  = alpha(y)*qeast;
      const double qlat_west  = alpha(y)*qwest;

      wtdnew(x,y) = wtd(x,y) + (qlat_north(x,y)+qlat_south(x,y)+qlat_east(x,y)+qlat_west(x,y))    //TODO: Check all of your signs! I'm not totally sure if it should be - or + here! 
      
      wtdnew(x,y+1) = wtd(x,y+1) - qlat_north;
      wtdnew(x,y-1) = wtd(x,y-1) - qlat_south;
      wtdnew(x-1,y) = wtd(x-1,y) - qlat_west;
      wtdnew(x+1,y) = wtd(x+1,y) - qlat_east;
    }



  }


}