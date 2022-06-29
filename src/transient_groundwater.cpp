#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <omp.h>
#include <array>
#include <experimental/source_location>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>
#include <petscts.h>

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

typedef enum { VAR_CONSERVATIVE, VAR_NONCONSERVATIVE, VAR_TRANSIENTVAR } VarMode;
static const char* const VarModes[] = {"CONSERVATIVE", "NONCONSERVATIVE", "TRANSIENTVAR", "VarMode", "VAR_", NULL};

double depthIntegratedTransmissivity(const double wtd_T, const double fdepth, const double ksat) {
  constexpr double shallow = 1.5;
  // Global soil datasets include information for shallow soils.
  // if the water table is deeper than this, the permeability
  // of the soil sees an exponential decay with depth.
  if (fdepth <= 0) {
    // If the fdepth is zero, there is no water transmission below the surface
    // soil layer.
    // If it is less than zero, it is incorrect -- but no water transmission
    // also seems an okay thing to do in this case.
    return 0;
  } else if (wtd_T < -shallow) {  // Equation S6 from the Fan paper
    return std::max(0.0, fdepth * ksat * std::exp((wtd_T + shallow) / fdepth));
  } else if (wtd_T > 0) {
    // If wtd_T is greater than 0, max out rate of groundwater movement
    // as though wtd_T were 0. The surface water will get to move in
    // FillSpillMerge.
    return std::max(0.0, ksat * (0 + shallow + fdepth));
  } else {                                                    // Equation S4 from the Fan paper
    return std::max(0.0, ksat * (wtd_T + shallow + fdepth));  // max because you can't have a negative transmissivity.
  }
}

static PetscErrorCode FormInitialSolution(DM da, Vec X, ArrayPack& arp) {
  PetscInt i, j, xs, ys, xm, ym, Mx, My;
  PetscScalar** x;
  DMDAGetInfo(
      da,
      PETSC_IGNORE,
      &Mx,
      &My,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE,
      PETSC_IGNORE);
  /*
     Get pointers to vector data
  */
  DMDAVecGetArray(da, X, &x);
  /*
     Get local grid boundaries
  */
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  /*
     Compute function over the locally owned part of the grid
  */
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      if (arp.land_mask(i, j) == 0) {
        x[i][j] = arp.topo(i, j) + 0.0;
      } else {
        x[i][j] = arp.topo(i, j) + arp.wtd(i, j);
      }
    }
  }
  /*
     Restore vectors
  */
  DMDAVecRestoreArray(da, X, &x);
  return 0;
}

static PetscErrorCode IFunction_TransientVar(DM da, PetscReal t, Vec X, Vec Cdot, Vec F, void* ctx) {
  auto* user_context = static_cast<AppCtx*>(ctx);
  da                 = user_context->da;
  // const PetscScalar* cdot;
  PetscInt xs, ys, xm, ym;
  Vec localX;
  Vec localC;
  Vec localF;

  // PetscScalar* f;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **my_h, **my_porosity,
      **my_topo, **x, **f, **my_rech, **my_area, **cdot;
  // const PetscScalar *x;
  std::cout << "doing the ifunction" << std::endl;

  PetscCall(DMDAVecGetArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecGetArray(da, user_context->S, &starting_storativity));
  PetscCall(DMDAVecGetArray(da, user_context->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecGetArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecGetArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecGetArray(da, user_context->porosity, &my_porosity));
  PetscCall(DMDAVecGetArray(da, user_context->h, &my_h));
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecGetArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecGetArray(da, user_context->cell_area, &my_area));

  DMGetLocalVector(da, &localX);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
  DMDAVecGetArrayRead(da, localX, &x);

  DMGetLocalVector(da, &localF);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, localF);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, localF);
  DMDAVecGetArrayRead(da, localF, &f);

  DMGetLocalVector(da, &localC);
  DMGlobalToLocalBegin(da, Cdot, INSERT_VALUES, localC);
  DMGlobalToLocalEnd(da, Cdot, INSERT_VALUES, localC);
  DMDAVecGetArrayRead(da, localC, &cdot);

  // VecGetArrayRead(X, &x);
  //  VecGetArrayRead(Cdot, &cdot);
  // PetscCall(DMDAVecGetArray(da, F, &f));

  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      //    std::cout<<"i "<<i<<" j "<<j<<" x "<<x[i][j]<<" cdot "<<cdot[i][j]<<" f "<<f[i][j]<<std::endl;
      const PetscScalar u = x[i][j];
      if (my_mask[i][j] == 0) {
        f[i][j] = u;
      } else {
        const PetscScalar ux_E = (x[i][j + 1] - x[i][j]);
        const PetscScalar ux_W = (x[i][j] - x[i][j - 1]);
        const PetscScalar uy_N = (x[i + 1][j] - x[i][j]);
        const PetscScalar uy_S = (x[i][j] - x[i - 1][j]);
        const PetscScalar e_E =
            2. / (1. / depthIntegratedTransmissivity(x[i][j] - my_topo[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(
                           x[i][j + 1] - my_topo[i][j + 1], my_fdepth[i][j + 1], my_ksat[i][j + 1]));
        const PetscScalar e_W =
            2. / (1. / depthIntegratedTransmissivity(x[i][j] - my_topo[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(
                           x[i][j - 1] - my_topo[i][j - 1], my_fdepth[i][j - 1], my_ksat[i][j - 1]));
        const PetscScalar e_N =
            2. / (1. / depthIntegratedTransmissivity(x[i][j] - my_topo[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(
                           x[i + 1][j] - my_topo[i + 1][j], my_fdepth[i + 1][j], my_ksat[i + 1][j]));
        const PetscScalar e_S =
            2. / (1. / depthIntegratedTransmissivity(x[i][j] - my_topo[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(
                           x[i - 1][j] - my_topo[i - 1][j], my_fdepth[i - 1][j], my_ksat[i - 1][j]));

        const PetscScalar uxx = (e_W * ux_W - e_E * ux_E) / (cellsize_ew[i][j] * cellsize_ew[i][j]);
        const PetscScalar uyy = (e_S * uy_S - e_N * uy_N) / (user_context->cellsize_NS * user_context->cellsize_NS);
        // TODO double check which cellsize is which

        starting_storativity[i][j] = updateEffectiveStorativity(
            my_h[i][j], x[i][j] - my_topo[i][j], my_porosity[i][j], starting_storativity[i][j]);

        double recharge = add_recharge(my_rech[i][j], my_h[i][j], my_porosity[i][j], my_area[i][j], 0);
        f[i][j] = x[i][j] + 10;  //(uxx + uyy) * (user_context->timestep / starting_storativity[i][j]) + u - recharge;
        //  std::cout<<"i "<<i<<" j "<<j<<" f "<<f[i][j]<<" x "<<x[i][j]<<std::endl;
      }
      //    std::cout<<"i "<<i<<" j "<<j<<" x "<<x[i][j]<<" cdot "<<cdot[i][j]<<" f "<<f[i][j]<<std::endl;
    }
  }
  std::cout << "done" << std::endl;
  std::cout << "x in the func " << x[10][10] << " f " << f[10][10] << std::endl;
  // f[0] = cdot[0] + PetscExpScalar(x[0]);
  // f[1] = cdot[1] - PetscExpScalar(x[0]);

  PetscCall(DMDAVecRestoreArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecRestoreArray(da, user_context->S, &starting_storativity));
  PetscCall(DMDAVecRestoreArray(da, user_context->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecRestoreArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecRestoreArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecRestoreArray(da, user_context->porosity, &my_porosity));
  PetscCall(DMDAVecRestoreArray(da, user_context->h, &my_h));
  PetscCall(DMDAVecRestoreArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecRestoreArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecRestoreArray(da, user_context->cell_area, &my_area));

  DMDAVecRestoreArrayRead(da, localX, &x);
  DMRestoreLocalVector(da, &localX);
  DMDAVecRestoreArrayRead(da, localF, &f);
  DMRestoreLocalVector(da, &localF);
  DMDAVecRestoreArrayRead(da, localC, &cdot);
  DMRestoreLocalVector(da, &localC);

  // VecRestoreArrayRead(Cdot, &cdot);
  std::cout << "return" << std::endl;
  return 0;
}

static PetscErrorCode TransientVar(TS ts, Vec X, Vec C, void* ctx) {
  std::cout << "doing the transient var" << std::endl;
  VecCopy(X, C);
  VecExp(C);
  return 0;
}

void PETSC_CHECK(
    const PetscErrorCode err,
    const std::experimental::source_location location = std::experimental::source_location::current()) {
  if (err) {
    throw std::runtime_error(
        "Petsc exception: " + std::to_string(err) + " at " + location.file_name() + ":" +
        std::to_string(location.line()));
  }
}

// get corners of arrays for individual processors
std::tuple<PetscInt, PetscInt, PetscInt, PetscInt> get_corners(const DM da) {
  PetscInt xs, ys, xm, ym;
  PETSC_CHECK(DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr));
  return {xs, ys, xm, ym};
}

extern PetscErrorCode MyTSMonitor(TS, PetscInt, PetscReal, Vec, void*);
extern PetscErrorCode MySNESMonitor(SNES, PetscInt, PetscReal, PetscViewerAndFormat*);

//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

void set_starting_values(Parameters& params, ArrayPack& arp) {
  // no pragma because we're editing arp.total_loss_to_ocean
  // check to see if there is any non-zero water table in ocean
  // cells, and if so, record these values as changes to the ocean.
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask(x, y) == 0.f) {
        arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y];
        arp.wtd(x, y) = 0.;
      } else {
        double rech_count = add_recharge(
            arp.rech(x, y),
            arp.wtd(x, y),
            arp.porosity(x, y),
            arp.cell_area[y],
            1);  // TODO: update how total_added_recharge gets counted
      }
    }
  }

#pragma omp parallel for default(none) shared(arp, params) collapse(2)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask(x, y) == 0.f || arp.wtd(x, y) > 0) {
        // in the ocean, we set several arrays to default values
        arp.effective_storativity(x, y) = 1.;
      } else {
        arp.effective_storativity(x, y) = arp.porosity(x, y);
      }
    }
  }
}

int update(Parameters& params, ArrayPack& arp, AppCtx& user_context, DMDA_Array_Pack& dmdapack) {
  TS ts;
  PetscInt xs, ys, xm, ym;
  PetscBool usemonitor = PETSC_TRUE;
  SNES ts_snes;
  PetscReal ftime;
  PetscInt steps; /* iterations for convergence */
  PetscViewerAndFormat* vf;
  PetscScalar sum;

  // compute any starting values needed for arrays
  set_starting_values(params, arp);

  std::cout << "line 137" << std::endl;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetDM(ts, user_context.da);
  TSSetProblemType(ts, TS_NONLINEAR);
  DMDAGetCorners(user_context.da, &xs, &ys, NULL, &xm, &ym, NULL);

  // values for storativity are reset each time; and wtd and recharge change from one timestep to the next, so set
  // these here
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.mask[i][j]        = arp.land_mask(i, j);
      dmdapack.S[i][j]           = arp.effective_storativity(i, j);
      dmdapack.cellsize_EW[i][j] = arp.cellsize_e_w_metres[j];
      dmdapack.fdepth_vec[i][j]  = arp.fdepth(i, j);
      dmdapack.ksat_vec[i][j]    = arp.ksat(i, j);
      dmdapack.porosity[i][j]    = arp.porosity(i, j);
      dmdapack.h[i][j]           = arp.wtd(i, j);
      dmdapack.topo_vec[i][j]    = arp.topo(i, j);
      dmdapack.rech_vec[i][j]    = arp.rech(i, j);
      dmdapack.cell_area[i][j]   = arp.cell_area[j];
    }
  }
  std::cout << "line 152 " << dmdapack.x[10][10] << std::endl;

  // DMDASNESSetFunctionLocal(user_context.da, INSERT_VALUES, (PetscErrorCode(*)(DMDALocalInfo*, void*, void*,
  // void*))FormFunctionLocal, &user_context);
  FormInitialSolution(user_context.da, user_context.x, arp);
  DMTSSetIFunctionLocal(user_context.da, IFunction_TransientVar, (void*)&user_context);
  // DMTSSetTransientVariable(user_context.da, TransientVar, (void*)&user_context);

  TSSetMaxTime(ts, 1.0);
  if (usemonitor) {
    TSMonitorSet(ts, MyTSMonitor, 0, 0);
  }
  std::cout << "line 158" << std::endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //  TSSetType(ts, TSBEULER);
  //  TSGetSNES(ts, &ts_snes);
  //  if (usemonitor) {
  //    PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
  //    SNESMonitorSet(
  //        ts_snes,
  //        (PetscErrorCode(*)(SNES, PetscInt, PetscReal, void*))MySNESMonitor,
  //        vf,
  //        (PetscErrorCode(*)(void**))PetscViewerAndFormatDestroy);
  //  }
  std::cout << "line 173" << std::endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSetTimeStep(ts, .1);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);
  //  TSSetSolution(ts, user_context.x);
  std::cout << "line 182" << std::endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSetFromOptions(ts);
  std::cout << "line 190" << std::endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSSolve(ts, user_context.x);
  std::cout << "line 196 " << dmdapack.x[10][10] << std::endl;

  TSGetSolveTime(ts, &ftime);
  std::cout << "line 199" << std::endl;

  TSGetStepNumber(ts, &steps);
  std::cout << "line 202" << std::endl;

  VecViewFromOptions(user_context.x, NULL, "-final_sol");
  std::cout << "line 205" << std::endl;

  // // VecCreateSeq(PETSC_COMM_SELF, 2, &U);
  // // VecSetValue(U, 0, 2., INSERT_VALUES);
  // // VecSetValue(U, 1, 1., INSERT_VALUES);
  // // std::cout << "line 187" << std::endl;
  // // VecLog(U);
  // // DMTSSetIFunction(user_context.da, IFunction_TransientVar, NULL);
  // // DMTSSetTransientVariable(user_context.da, TransientVar, NULL);
  // // std::cout << "line 191" << std::endl;
  // // TSSetMaxTime(ts, 1.);
  // // TSSetFromOptions(ts);
  // // TSSolve(ts, U);
  // // std::cout << "line 195" << std::endl;
  // // VecExp(U);
  //
  // // VecView(U, PETSC_VIEWER_STDOUT_WORLD);
  // // VecSum(U, &sum);
  // // PetscPrintf(PETSC_COMM_WORLD, "Conservation error %g\n", (double)PetscRealPart(sum - 3.));
  // VecDestroy(&U);
  std::cout << "line 209" << std::endl;

  // copy the result back into the wtd array
  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      arp.wtd(i, j) = dmdapack.x[i][j] - arp.topo(i, j);
      if (arp.land_mask(i, j) == 0.f) {
        arp.total_loss_to_ocean +=
            arp.wtd(i, j) * arp.cell_area[j];  // could it be that because ocean cells are just set = x in the formula,
                                               // that loss to/gain from ocean is not properly recorded?
        arp.wtd(i, j) = 0.;
      }
    }
  }

  TSDestroy(&ts);
  return 0;
}

/*
   FormRHS - Forms constant RHS for the problem.

   Input Parameters:
   user - user-defined application context
   B - RHS vector

   Output Parameter:
   B - vector
 */
// static PetscErrorCode FormRHS(AppCtx* user_context, DM da, Vec B, ArrayPack& arp) {
//   PetscInt Mx, My;
//   PetscScalar** b;
//
//   DMDAGetInfo(
//       da,
//       PETSC_IGNORE,
//       &Mx,
//       &My,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE,
//       PETSC_IGNORE);
//
//   DMDAVecGetArray(da, B, &b);
//
//   const auto [xs, ys, xm, ym] = get_corners(da);
//   for (auto j = ys; j < ys + ym; j++) {
//     for (auto i = xs; i < xs + xm; i++) {
//       if (arp.land_mask(i, j) == 0) {
//         b[i][j] = arp.topo(i, j) + 0.0;
//       } else {
//         b[i][j] = arp.topo(i, j) + arp.wtd(i, j) +
//                   add_recharge(arp.rech(i, j), arp.wtd(i, j), arp.porosity(i, j), arp.cell_area[j], 0, arp);
//       }
//     }
//   }
//   DMDAVecRestoreArray(da, B, &b);
//   return 0;
// }
//
///* ------------------------------------------------------------------- */

PetscErrorCode MyTSMonitor(TS ts, PetscInt step, PetscReal ptime, Vec v, void* ctx) {
  PetscReal norm;
  MPI_Comm comm;
  VecNorm(v, NORM_2, &norm);
  PetscObjectGetComm((PetscObject)ts, &comm);
  if (step > -1) { /* -1 is used to indicate an interpolated value */
    PetscPrintf(comm, "timestep %" PetscInt_FMT " time %g norm %g\n", step, (double)ptime, (double)norm);
  }
  return 0;
}
/*
   MySNESMonitor - illustrate how to set user-defined monitoring routine for SNES.
   Input Parameters:
     snes - the SNES context
     its - iteration number
     fnorm - 2-norm function value (may be estimated)
     ctx - optional user-defined context for private data for the
         monitor routine, as set by SNESMonitorSet()
 */
PetscErrorCode MySNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, PetscViewerAndFormat* vf) {
  SNESMonitorDefaultShort(snes, its, fnorm, vf);
  return 0;
}
/*TEST
    test:
      args: -da_grid_x 20 -ts_max_time 3 -ts_dt 1e-1 -ts_theta_initial_guess_extrapolate 0 -ts_monitor
-ksp_monitor_short requires: !single test: suffix: 2 args: -da_grid_x 20 -ts_max_time 0.11 -ts_dt 1e-1 -ts_type glle
-ts_monitor -ksp_monitor_short requires: !single test: suffix: glvis_da_2d_vect args: -usemonitor 0 -da_grid_x 20
-ts_max_time 0.3 -ts_dt 1e-1 -ts_type glle -final_sol glvis: -viewer_glvis_dm_da_bs 2,0 requires: !single test: suffix:
glvis_da_2d_vect_ll args: -usemonitor 0 -da_grid_x 20 -ts_max_time 0.3 -ts_dt 1e-1 -ts_type glle -final_sol glvis:
-viewer_glvis_dm_da_bs 2,0 -viewer_glvis_dm_da_ll requires: !single TEST*/

}  // namespace FanDarcyGroundwater
