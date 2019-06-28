#include "ArrayPack.hpp"
#include <cassert>

void ArrayPack::check() const {
  assert( topo_start.width()==ksat.         width() && topo_start.height()==ksat.         height() );
  assert( topo_start.width()==land.         width() && topo_start.height()==land.         height() );
  assert( topo_start.width()==wtd.          width() && topo_start.height()==wtd.          height() );
  assert( topo_start.width()==topo_start.   width() && topo_start.height()==topo_start.   height() );
 // assert( topo_start.width()==rech_start.   width() && topo_start.height()==rech_start.   height() );
 // assert( topo_start.width()==fslope_start. width() && topo_start.height()==fslope_start. height() );
 // assert( topo_start.width()==temp_start.   width() && topo_start.height()==temp_start.   height() );
//  assert( topo_start.width()==topo.         width() && topo_start.height()==topo.         height() );
  assert( topo_start.width()==rech.         width() && topo_start.height()==rech.         height() );
  assert( topo_start.width()==fdepth.       width() && topo_start.height()==fdepth.       height() );
  assert( topo_start.width()==temp.         width() && topo_start.height()==temp.         height() );
  assert( topo_start.width()==done_old.     width() && topo_start.height()==done_old.     height() );
  assert( topo_start.width()==done_new.     width() && topo_start.height()==done_new.     height() );
  assert( topo_start.width()==head.         width() && topo_start.height()==head.         height() );
  assert( topo_start.width()==kcell.        width() && topo_start.height()==kcell.        height() );


  if(fslope_end.size()>0){
    assert( topo_start.width()==fslope_end.width()   && topo_start.height()==fslope_end.height()   );
    assert( topo_start.width()==topo_end.  width()   && topo_start.height()==topo_end.  height()   );
    assert( topo_start.width()==temp_end.  width()   && topo_start.height()==temp_end.  height()   );
    assert( topo_start.width()==rech_end.  width()   && topo_start.height()==rech_end.  height()   );
  }
}