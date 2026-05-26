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
      dmdapack.cellsize_EW_squared[j][i] = arp.cellsize_e_w_metres[j] * arp.cellsize_e_w_metres[j];
      dmdapack.mask[j][i]                = arp.land_mask(i, j);
      dmdapack.porosity_vec[j][i]        = arp.porosity(i, j);
    }
  }
}

// Populate the global topo/fdepth/ksat vecs from arp and scatter to local ghost vectors.
// Must be called while these global vecs are NOT under DMDAVecGetArray (i.e., before or after
// DMDA_Array_Pack holds them — which it does NOT, by design).
void scatter_static_fields(AppCtx& user_context, ArrayPack& arp) {
  const auto [xs, ys, xm, ym] = get_corners(user_context.da);
  PetscScalar **topo_arr, **fdepth_arr, **ksat_arr;

  DMDAVecGetArray(user_context.da, user_context.topo_vec, &topo_arr);
  DMDAVecGetArray(user_context.da, user_context.fdepth_vec, &fdepth_arr);
  DMDAVecGetArray(user_context.da, user_context.ksat_vec, &ksat_arr);
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      topo_arr[j][i]   = arp.topo(i, j);
      fdepth_arr[j][i] = arp.fdepth(i, j);
      ksat_arr[j][i]   = arp.ksat(i, j);
    }
  }
  DMDAVecRestoreArray(user_context.da, user_context.topo_vec, &topo_arr);
  DMDAVecRestoreArray(user_context.da, user_context.fdepth_vec, &fdepth_arr);
  DMDAVecRestoreArray(user_context.da, user_context.ksat_vec, &ksat_arr);

  // The DMDA's internal PetscSF is shared across all GlobalToLocal operations on the same DM.
  // Overlapping Begin calls (Begin A, Begin B, End A, End B) confuse the SF state machine;
  // each pair must be completed sequentially.
  DMGlobalToLocalBegin(user_context.da, user_context.topo_vec,   INSERT_VALUES, user_context.topo_local);
  DMGlobalToLocalEnd(user_context.da, user_context.topo_vec,     INSERT_VALUES, user_context.topo_local);
  DMGlobalToLocalBegin(user_context.da, user_context.fdepth_vec, INSERT_VALUES, user_context.fdepth_local);
  DMGlobalToLocalEnd(user_context.da, user_context.fdepth_vec,   INSERT_VALUES, user_context.fdepth_local);
  DMGlobalToLocalBegin(user_context.da, user_context.ksat_vec,   INSERT_VALUES, user_context.ksat_local);
  DMGlobalToLocalEnd(user_context.da, user_context.ksat_vec,     INSERT_VALUES, user_context.ksat_local);
}
