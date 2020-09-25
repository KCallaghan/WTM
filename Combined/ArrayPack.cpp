#include "ArrayPack.hpp"
#include <cassert>

void ArrayPack::check() const {
  assert( topo.width()==ksat.         width() &&     topo.height()==ksat.         height() );
  assert( topo.width()==land_mask.    width() &&     topo.height()==land_mask.    height() );
  assert( topo.width()==wtd.          width() &&     topo.height()==wtd.          height() );
  assert( topo.width()==fdepth.       width() &&     topo.height()==fdepth.       height() );
  assert( topo.width()==precip.       width() &&     topo.height()==precip.       height() );
  assert( topo.width()==fdepth.       width() &&     topo.height()==fdepth.       height() );
  assert( topo.width()==starting_evap.width() &&     topo.height()==starting_evap.height() );
  assert( topo.width()==head.         width() &&     topo.height()==head.         height() );
  assert( topo.width()==open_water_evap. width() &&     topo.height()==open_water_evap. height() );
  assert( topo.width()==runoff.       width() &&     topo.height()==runoff.       height() );

  if(fdepth_end.size()>0){
    assert( topo.width()==fdepth_end.width()  &&       topo.height()==fdepth_end.   height()   );
    assert( topo.width()==topo_end.  width()  &&       topo.height()==topo_end.     height()   );
    assert( topo.width()==temp_end.  width()  &&       topo.height()==temp_end.     height()   );
    assert( topo.width()==precip_end.width()  &&       topo.height()==precip_end.   height()   );
  }
}