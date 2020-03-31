#include "ArrayPack.hpp"
#include <cassert>

void ArrayPack::check() const {
  assert( topo.width()==ksat.         width() && topo.height()==ksat.         height() );
  assert( topo.width()==land_mask.    width() && topo.height()==land_mask.    height() );
  assert( topo.width()==wtd.          width() && topo.height()==wtd.          height() );
  assert( topo.width()==fdepth.       width() && topo.height()==fdepth.       height() );
  assert( topo.width()==temp.         width() && topo.height()==temp.         height() );
  assert( topo.width()==slope.        width() && topo.height()==slope.        height() );
  assert( topo.width()==precip.       width() && topo.height()==precip.       height() );
  assert( topo.width()==fdepth.       width() && topo.height()==fdepth.       height() );
  assert( topo.width()==temp.         width() && topo.height()==temp.         height() );
  assert( topo.width()==starting_evap.width() && topo.height()==starting_evap.height() );
  assert( topo.width()==relhum.       width() && topo.height()==relhum.       height() );
  assert( topo.width()==head.         width() && topo.height()==head.         height() );
  assert( topo.width()==kcell.        width() && topo.height()==kcell.        height() );
  assert( topo.width()==evap.         width() && topo.height()==evap.         height() );
  assert( topo.width()==e_sat.        width() && topo.height()==e_sat.        height() );
  assert( topo.width()==e_a.          width() && topo.height()==e_a.          height() );
  assert( topo.width()==surface_evap. width() && topo.height()==surface_evap. height() );
  assert( topo.width()==runoff.       width() && topo.height()==runoff.       height() );

  if(fdepth_end.size()>0){
    assert( topo.width()==fdepth_end.width()  && topo.height()==fdepth_end.   height()   );
    assert( topo.width()==topo_end.  width()  && topo.height()==topo_end.     height()   );
    assert( topo.width()==temp_end.  width()  && topo.height()==temp_end.     height()   );
    assert( topo.width()==precip_end.width()  && topo.height()==precip_end.   height()   );
  }
}