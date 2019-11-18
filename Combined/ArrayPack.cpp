#include "ArrayPack.hpp"
#include <cassert>

void ArrayPack::check() const {
  assert( topo.width()==ksat.         width() && topo.height()==ksat.         height() );
  assert( topo.width()==land_mask.         width() && topo.height()==land_mask.         height() );
  assert( topo.width()==wtd.          width() && topo.height()==wtd.          height() );
//  assert( topo.width()==topo.   width() && topo.height()==topo.   height() );
//  assert( topo.width()==rech_start.   width() && topo.height()==rech_start.   height() );
//  assert( topo.width()==fslope_start. width() && topo.height()==fslope_start. height() );
//  assert( topo.width()==temp_start.   width() && topo.height()==temp_start.   height() );
//  assert( topo.width()==topo.         width() && topo.height()==topo.         height() );
  assert( topo.width()==precip.         width() && topo.height()==precip.         height() );
  assert( topo.width()==fdepth.       width() && topo.height()==fdepth.       height() );
  assert( topo.width()==temp.         width() && topo.height()==temp.         height() );
 // assert( topo.width()==done_old.     width() && topo.height()==done_old.     height() );
 // assert( topo.width()==done_new.     width() && topo.height()==done_new.     height() );
 // assert( topo.width()==head.         width() && topo.height()==head.         height() );
 // assert( topo.width()==kcell.        width() && topo.height()==kcell.        height() );


  if(fslope_end.size()>0){
    assert( topo.width()==fslope_end.width()   && topo.height()==fslope_end.height()   );
    assert( topo.width()==topo_end.  width()   && topo.height()==topo_end.  height()   );
    assert( topo.width()==temp_end.  width()   && topo.height()==temp_end.  height()   );
    assert( topo.width()==rech_end.  width()   && topo.height()==rech_end.  height()   );
  }
}