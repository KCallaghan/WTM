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
      dmdapack.cellsize_EW[i][j] = arp.cellsize_e_w_metres[j];
      dmdapack.mask[i][j]        = arp.land_mask(i, j);
      dmdapack.fdepth_vec[i][j]  = arp.fdepth(i, j);
      dmdapack.ksat_vec[i][j]    = arp.ksat(i, j);
      dmdapack.porosity[i][j]    = arp.porosity(i, j);
      dmdapack.topo_vec[i][j]    = arp.topo(i, j);
    }
  }
}
