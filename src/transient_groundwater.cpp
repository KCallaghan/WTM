#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <omp.h>
#include <array>
#include <chrono>
#include <experimental/source_location>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

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

// declare functions
static PetscErrorCode FormRHS(AppCtx*, DM, Vec);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec);
static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);

//////////////////////
// PUBLIC FUNCTIONS //
//////////////////////

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

void set_starting_values(Parameters& params, ArrayPack& arp) {
  // no pragma because we're editing arp.total_loss_to_ocean
  // check to see if there is any non-zero water table in ocean
  // cells, and if so, record these values as changes to the ocean.
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask(x, y) == 0.f) {
        if (arp.wtd(x, y) > 0)
          arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y];
        else
          arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y] * arp.porosity(x, y);
        arp.wtd(x, y) = 0.;
      } else {
        double rech_count = arp.rech(x, y);
        if (arp.wtd(x, y) >= 0 && arp.wtd(x, y) + arp.rech(x, y) < 0)
          rech_count = -arp.wtd(x, y);

        arp.total_added_recharge += rech_count * arp.cell_area[y];
      }
    }
  }
}

int update(Parameters& params, ArrayPack& arp, AppCtx& user_context, DMDA_Array_Pack& dmdapack) {
  PetscInt its;                // iterations for convergence
  SNESConvergedReason reason;  // Check convergence

  // compute any starting values needed for arrays
  set_starting_values(params, arp);

  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user_context.da);

//  values for storativity are reset each time; and recharge changes from one timestep to the next, so set these here
#pragma omp parallel for default(none) shared(arp, ys, ym, xs, xm, dmdapack, params) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.rech_vec[j][i]     = add_recharge(arp.rech(i, j), arp.wtd(i, j), arp.porosity(i, j));
      dmdapack.head[j][i]         = arp.wtd(i, j) + arp.topo(i, j);
      dmdapack.starting_wtd[j][i] = arp.wtd(i, j);
    }
  }

  // Set local function evaluation routine
  DMDASNESSetFunctionLocal(
      user_context.da,
      INSERT_VALUES,
      (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal,
      &user_context);

  // Evaluate initial guess
  FormInitialGuess(&user_context, user_context.da, user_context.x);

  // set the RHS
  FormRHS(&user_context, user_context.da, user_context.b);
  // Solve nonlinear system
  SNESSolve(user_context.snes, user_context.b, user_context.x);

  SNESGetIterationNumber(user_context.snes, &its);
  SNESGetConvergedReason(user_context.snes, &reason);

  PetscPrintf(
      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);

  // copy the result back into the wtd array
  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      arp.wtd(i, j) = dmdapack.x[j][i] - arp.topo(i, j);
      if (arp.land_mask(i, j) == 0.f) {
        if (arp.wtd(i, j) > 0)
          arp.total_loss_to_ocean +=
              arp.wtd(i, j) * arp.cell_area[j];  // could it be that because ocean cells are just set = x in the
                                                 // formula, that loss to/gain from ocean is not properly recorded?
        else
          arp.total_loss_to_ocean += arp.wtd(i, j) * arp.cell_area[j] * arp.porosity(i, j);
        arp.wtd(i, j) = 0.;
      }
    }
  }

  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormInitialGuess - Forms initial approximation.

   Input Parameters:
   user - user-defined application context
   X - vector

   Output Parameter:
   X - vector
 */
static PetscErrorCode FormInitialGuess(AppCtx* user_context, DM da, Vec X) {
  PetscScalar **x, **my_head;

  DMDAVecGetArray(da, X, &x);
  DMDAVecGetArray(da, user_context->head, &my_head);

  const auto [xs, ys, xm, ym] = get_corners(da);

#pragma omp parallel for default(none) shared(my_head, ys, ym, xs, xm, x) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      x[j][i] = my_head[j][i];  // when land mask == 0, both topo and wtd have already been set to 0
                                // elsewhere, so no need for another if statement here
    }
  }

  DMDAVecRestoreArray(da, X, &x);
  DMDAVecRestoreArray(da, user_context->head, &my_head);
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
static PetscErrorCode FormRHS(AppCtx* user_context, DM da, Vec B) {
  PetscScalar **b, **my_head;

  DMDAVecGetArray(da, B, &b);
  DMDAVecGetArray(da, user_context->head, &my_head);

  const auto [xs, ys, xm, ym] = get_corners(da);

#pragma omp parallel for default(none) shared(ys, ym, xs, xm, b, my_head) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      b[j][i] = my_head[j][i];  // when land mask == 0, both topo and wtd have already been set to 0
                                // elsewhere, so no need for another if statement here
    }
  }
  DMDAVecRestoreArray(da, B, &b);
  DMDAVecRestoreArray(da, user_context->head, &my_head);

  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x).
 */
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user_context) {
  DM da = user_context->da;
  PetscScalar **cellsize_ew_sq, **my_mask, **my_fdepth, **my_ksat, **my_topo, **my_rech, **my_T, **my_starting_wtd,
      **my_porosity;

  /*
    Compute function over the locally owned part of the grid
 */
  PetscCall(DMDAVecGetArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecGetArray(da, user_context->cellsize_EW_squared, &cellsize_ew_sq));
  PetscCall(DMDAVecGetArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecGetArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecGetArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecGetArray(da, user_context->T_vec, &my_T));
  PetscCall(DMDAVecGetArray(da, user_context->porosity_vec, &my_porosity));
  PetscCall(DMDAVecGetArray(da, user_context->starting_wtd, &my_starting_wtd));

#pragma omp parallel for default(none) shared(info, my_T, x, my_topo, my_fdepth, my_ksat) collapse(2)
  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      my_T[j][i] = 1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]);
    }
  }

#pragma omp parallel for default(none)                                                                              \
    shared(info, cellsize_ew_sq, x, my_T, my_mask, my_rech, user_context, my_porosity, my_starting_wtd, my_topo, f) \
        collapse(2)
  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      if (my_mask[j][i] == 0) {
        f[j][i] = 0;
      } else {
        double this_x          = x[j][i];
        double this_T          = my_T[j][i];
        const PetscScalar ux_E = (x[j][i + 1] - this_x);
        const PetscScalar ux_W = (this_x - x[j][i - 1]);
        const PetscScalar uy_N = (x[j + 1][i] - this_x);
        const PetscScalar uy_S = (this_x - x[j - 1][i]);
        const PetscScalar e_E  = 2. / (this_T + my_T[j][i + 1]);
        const PetscScalar e_W  = 2. / (this_T + my_T[j][i - 1]);
        const PetscScalar e_N  = 2. / (this_T + my_T[j + 1][i]);
        const PetscScalar e_S  = 2. / (this_T + my_T[j - 1][i]);

        const PetscScalar uxx = (e_W * ux_W - e_E * ux_E) / user_context->cellsize_NS_squared;
        const PetscScalar uyy = (e_S * uy_S - e_N * uy_N) / cellsize_ew_sq[j][i];

        double my_storativity =
            updateEffectiveStorativity(my_starting_wtd[j][i], this_x - my_topo[j][i], my_porosity[j][i]);

        f[j][i] = (uxx + uyy) * user_context->deltat / my_storativity + this_x - my_rech[j][i];
        // my_rech is converted to appropriate recharge for this timestep and starting water
        // table outside of the solve.
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecRestoreArray(da, user_context->cellsize_EW_squared, &cellsize_ew_sq));
  PetscCall(DMDAVecRestoreArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecRestoreArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecRestoreArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecRestoreArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecRestoreArray(da, user_context->T_vec, &my_T));
  PetscCall(DMDAVecRestoreArray(da, user_context->porosity_vec, &my_porosity));
  PetscCall(DMDAVecRestoreArray(da, user_context->starting_wtd, &my_starting_wtd));

  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

}  // namespace FanDarcyGroundwater
