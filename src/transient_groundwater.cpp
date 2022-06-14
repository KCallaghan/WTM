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

std::tuple<PetscInt, PetscInt, PetscInt, PetscInt> get_corners(const DM da) {
  PetscInt xs, ys, xm, ym;
  PETSC_CHECK(DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr));
  return {xs, ys, xm, ym};
}

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal().
*/

// User-defined routines

static PetscErrorCode FormRHS(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);

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
        double rech_count = add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);
      }
    }
  }

  //#pragma omp parallel for default(none) shared(arp, params) collapse(2)
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
  PetscInt its;               /* iterations for convergence */
  SNESConvergedReason reason; /* Check convergence */

  set_starting_values(params, arp);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DM; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user_context.da);

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.S[i][j]        = arp.effective_storativity(i, j);
      dmdapack.h[i][j]        = arp.wtd(i, j);
      dmdapack.rech_vec[i][j] = arp.rech(i, j);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set local function evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDASNESSetFunctionLocal(
      user_context.da,
      INSERT_VALUES,
      (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal,
      &user_context);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */
  FormInitialGuess(&user_context, user_context.da, user_context.x, arp);
  FormRHS(&user_context, user_context.da, user_context.b, arp);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESSolve(user_context.snes, user_context.b, user_context.x);
  SNESGetIterationNumber(user_context.snes, &its);
  SNESGetConvergedReason(user_context.snes, &reason);

  PetscPrintf(
      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.S[i][j] = updateEffectiveStorativity(
          arp.wtd(i, j), dmdapack.x[i][j] - arp.topo(i, j), arp.porosity(i, j), arp.effective_storativity(i, j));
    }
  }

  SNESSolve(user_context.snes, user_context.b, user_context.x);
  SNESGetIterationNumber(user_context.snes, &its);
  SNESGetConvergedReason(user_context.snes, &reason);

  PetscPrintf(
      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
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
static PetscErrorCode FormInitialGuess(AppCtx* user_context, DM da, Vec X, ArrayPack& arp) {
  PetscInt Mx, My;
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

  DMDAVecGetArray(da, X, &x);

  const auto [xs, ys, xm, ym] = get_corners(da);

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      if (arp.land_mask(i, j) == 0) {
        /* boundary conditions are all zero Dirichlet */
        x[i][j] = arp.topo(i, j) + 0.0;
      } else {
        x[i][j] = arp.topo(i, j) + arp.wtd(i, j);
      }
    }
  }

  DMDAVecRestoreArray(da, X, &x);
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
static PetscErrorCode FormRHS(AppCtx* user_context, DM da, Vec B, ArrayPack& arp) {
  PetscInt Mx, My;
  PetscScalar** b;

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

  DMDAVecGetArray(da, B, &b);

  const auto [xs, ys, xm, ym] = get_corners(da);
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      if (arp.land_mask(i, j) == 0) {
        b[i][j] = arp.topo(i, j) + 0.0;
      } else {
        b[i][j] = arp.topo(i, j) + arp.wtd(i, j) +
                  add_recharge(arp.rech(i, j), arp.wtd(i, j), arp.porosity(i, j), arp.cell_area[j], 0, arp);
      }
    }
  }
  DMDAVecRestoreArray(da, B, &b);
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x).
 */
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user_context) {
  DM da = user_context->da;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **my_h, **my_porosity,
      **my_topo, **my_rech, **my_cell_area;
  PetscInt Mx, My;

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
    Compute function over the locally owned part of the grid
 */
  PetscCall(DMDAVecGetArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecGetArray(da, user_context->S, &starting_storativity));
  PetscCall(DMDAVecGetArray(da, user_context->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecGetArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecGetArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecGetArray(da, user_context->porosity, &my_porosity));
  PetscCall(DMDAVecGetArray(da, user_context->h, &my_h));
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));

  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
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
        // std::cout<<"j "<<j<<" i "<<i<<" wtd "<<x[i][j]- my_topo[i][j]<<" fdepth "<<my_fdepth[i][j]<<" ksat
        // "<<my_ksat[i][j]<<" T "<<depthIntegratedTransmissivity(x[i][j]-
        // my_topo[i][j],my_fdepth[i][j],my_ksat[i][j])<<" e_E "<<e_E<<std::endl;
        //  starting_storativity[i][j] = updateEffectiveStorativity(
        //    my_h[i][j], x[i][j] - my_topo[i][j], my_porosity[i][j], starting_storativity[i][j]);

        f[i][j] = (uxx + uyy) * (user_context->timestep / starting_storativity[i][j]) + u;
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

  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

}  // namespace FanDarcyGroundwater
