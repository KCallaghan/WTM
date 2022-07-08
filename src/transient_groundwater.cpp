#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscts.h>

namespace FanDarcyGroundwater {

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

/* ------------------------------------------------------------------- */

/*
  FormInitialSolution - Forms initial approximation.

  Input Parameters:
  user - user-defined application context
  X - vector

  Output Parameter:
  X - vector
*/
static PetscErrorCode FormInitialSolution(DM da, Vec X, ArrayPack& arp) {
  PetscInt i, j, xs, ys, xm, ym;
  PetscScalar** x;

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
        x[j][i] = arp.topo(i, j) + 0.0;
      } else {
        x[j][i] = arp.topo(i, j) + arp.wtd(i, j);
      }
    }
  }
  std::cout << "x in the initial solution " << x[10][10] << std::endl;
  /*
     Restore vectors
  */
  DMDAVecRestoreArray(da, X, &x);
  return 0;
}

static PetscErrorCode FormIFunctionLocal(
    DMDALocalInfo* info,
    PetscReal ptime,
    PetscScalar** x,
    PetscScalar** Thickness,
    PetscScalar** f,
    void* ctx) {
  AppCtx* user_context = (AppCtx*)ctx;
  DM da                = user_context->da;
  PetscInt xs, ys, xm, ym;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **my_h, **my_porosity,
      **my_topo, **my_rech, **my_area;
  std::cout << "x in the func before " << x[10][10] << " f " << f[10][10] << " xdot " << xdot[10][10] << std::endl;

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

  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      const PetscScalar u = x[j][i];
      if (my_mask[j][i] == 0) {
        f[j][i] = u;
      } else {
        const PetscScalar ux_E = (x[j][i + 1] - x[j][i]);
        const PetscScalar ux_W = (x[j][i] - x[j][i - 1]);
        const PetscScalar uy_N = (x[j + 1][i] - x[j][i]);
        const PetscScalar uy_S = (x[j][i] - x[j - 1][i]);
        const PetscScalar e_E =
            2. / (1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]) +
                  1. / depthIntegratedTransmissivity(
                           x[j][i + 1] - my_topo[j][i + 1], my_fdepth[j][i + 1], my_ksat[j][i + 1]));
        const PetscScalar e_W =
            2. / (1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]) +
                  1. / depthIntegratedTransmissivity(
                           x[j][i - 1] - my_topo[j][i - 1], my_fdepth[j][i - 1], my_ksat[j][i - 1]));
        const PetscScalar e_N =
            2. / (1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]) +
                  1. / depthIntegratedTransmissivity(
                           x[j + 1][i] - my_topo[j + 1][i], my_fdepth[j + 1][i], my_ksat[j + 1][i]));
        const PetscScalar e_S =
            2. / (1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]) +
                  1. / depthIntegratedTransmissivity(
                           x[j - 1][i] - my_topo[j - 1][i], my_fdepth[j - 1][i], my_ksat[j - 1][i]));

        const PetscScalar uxx = (e_E * ux_E - e_W * ux_W) / (cellsize_ew[j][i] * cellsize_ew[j][i]);
        const PetscScalar uyy = (e_N * uy_N - e_S * uy_S) / (user_context->cellsize_NS * user_context->cellsize_NS);
        // TODO double check which cellsize is which

        starting_storativity[j][i] = updateEffectiveStorativity(
            my_h[j][i], x[j][i] - my_topo[j][i], my_porosity[j][i], starting_storativity[j][i]);

        double recharge = add_recharge(my_rech[j][i], my_h[j][i], my_porosity[j][i], my_area[j][i], 0);

        f[j][i] = Thickness[j][i] - (uxx + uyy + recharge);  //(uxx + uyy + recharge);
      }
    }
  }

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

  return 0;
}

static PetscErrorCode TransientVar(TS ts, Vec Head, Vec Thickness, void* ctx) {
  AppCtx* user_context = (AppCtx*)ctx;
  DM da                = user_context->da;
  // const PetscScalar* cdot;
  PetscInt xs, ys, xm, ym;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **H, **my_porosity, **my_topo,
      **my_rech, **my_area;
  const PetscScalar** h;
  PetscCall(DMDAVecGetArrayRead(da, Head, &h));
  PetscCall(DMDAVecGetArrayWrite(da, Thickness, &H));
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecGetArray(da, user_context->S, &starting_storativity));

  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      H[j][i] = (h[j][i] - my_topo[j][i]) * starting_storativity[j][i];
    }
  }

  PetscCall(DMDAVecRestoreArrayRead(da, Head, &h));
  PetscCall(DMDAVecRestoreArrayWrite(da, Thickness, &H));
  PetscCall(DMDAVecRestoreArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecRestoreArray(da, user_context->S, &starting_storativity));

  return 0;
}

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
  PetscInt steps;
  PetscReal ftime;
  DMDALocalInfo* info;
  PetscInt xs, ys, xm, ym;
  TSConvergedReason reason;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **my_h, **my_porosity,
      **my_topo, **my_rech, **my_area, **my_x, **my_f;

  set_starting_values(params, arp);

  TSSetProblemType(user_context.ts, TS_NONLINEAR);
  TSSetTime(user_context.ts, 0.0);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create user context, set problem data, create vector data structures.
   Also, compute the initial guess.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDAGetCorners(user_context.da, &xs, &ys, NULL, &xm, &ym, NULL);

  // values for storativity are reset each time; and wtd and recharge change from one timestep to the next, so set
  // these here

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.S[j][i]        = arp.effective_storativity(i, j);
      dmdapack.h[j][i]        = arp.wtd(i, j);
      dmdapack.rech_vec[j][i] = arp.rech(i, j);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create time integration context
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMSetApplicationContext(user_context.da, &user_context);
  DMDATSSetIFunctionLocal(
      user_context.da,
      INSERT_VALUES,
      (DMDATSIFunctionLocal)FormIFunctionLocal /*(info,ftime,dmdapack.x,dmdapack.xdot,dmdapack.f,&user_context)*/,
      (void*)&user_context);
  DMTSSetTransientVariable(user_context.da, TransientVar, (void*)&user_context);
  TSSetMaxSteps(user_context.ts, 10000);
  TSSetMaxTime(user_context.ts, user_context.maxtime);
  TSSetExactFinalTime(user_context.ts, TS_EXACTFINALTIME_STEPOVER);
  TSSetTimeStep(user_context.ts, user_context.timestep);
  TSSetFromOptions(user_context.ts);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  FormInitialSolution(user_context.da, user_context.x, arp);
  TSSetSolution(user_context.ts, user_context.x);
  TSSolve(user_context.ts, user_context.x);
  TSGetSolveTime(user_context.ts, &ftime);
  TSGetStepNumber(user_context.ts, &steps);
  TSGetConvergedReason(user_context.ts, &reason);

  PetscPrintf(
      PETSC_COMM_WORLD,
      "%s at time %g after %" PetscInt_FMT " steps\n",
      TSConvergedReasons[reason],
      (double)ftime,
      steps);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // copy the result back into the wtd array
  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      arp.wtd(i, j) = dmdapack.x[j][i] - arp.topo(i, j);
      if (arp.land_mask(i, j) == 0.f) {
        arp.total_loss_to_ocean +=
            arp.wtd(i, j) * arp.cell_area[j];  // could it be that because ocean cells are just set = x in the formula,
                                               // that loss to/gain from ocean is not properly recorded?
        arp.wtd(i, j) = 0.;
      }
    }
  }

  return 0;
}

}  // namespace FanDarcyGroundwater
