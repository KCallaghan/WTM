std::tuple<PetscInt, PetscInt, PetscInt, PetscInt> get_corners(const DM da) {
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr);
  return {xs, ys, xm, ym};
}

void populate_DMDA_array_pack(AppCtx& user_context, ArrayPack& arp, DMDA_Array_Pack& dmdapack) {
  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user_context.da);

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.cellsize_EW[j][i] = arp.cellsize_e_w_metres[j];
      dmdapack.mask[j][i]        = arp.land_mask(i, j);
      dmdapack.fdepth_vec[j][i]  = arp.fdepth(i, j);
      dmdapack.ksat_vec[j][i]    = arp.ksat(i, j);
      dmdapack.porosity[j][i]    = arp.porosity(i, j);
      dmdapack.topo_vec[j][i]    = arp.topo(i, j);
      //     std::cout<<" i "<<i<<" j "<<j<<" mask arp "<<arp.land_mask(i,j)<<" dmdapack
      //     "<<dmdapack.mask[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<" fdepth arp "<<arp.fdepth(i,j)<<"
      //     dmdapack "<<dmdapack.fdepth_vec[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<" ksat arp
      //     "<<arp.ksat(i,j)<<" dmdapack "<<dmdapack.ksat_vec[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<"
      //     porosity arp "<<arp.porosity(i,j)<<" dmdapack "<<dmdapack.porosity[i][j]<<std::endl;
    }
  }

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      //   std::cout<<" i "<<i<<" j "<<j<<" at creation mask "    <<arp.land_mask(i,j)<<"
      //   "<<dmdapack.mask[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<" at creation fdepth "  <<arp.fdepth(i,j)
      //   <<"    "<<dmdapack.fdepth_vec[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<" at creation ksat "
      //   <<arp.ksat(i,j)     <<"    "<<dmdapack.ksat_vec[i][j]<<std::endl; std::cout<<" i "<<i<<" j "<<j<<" at
      //   creation porosity "<<arp.porosity(i,j) <<"    "<<dmdapack.porosity[i][j]<<std::endl;
    }
  }
}
