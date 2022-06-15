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
static PetscErrorCode IFunction_Conservative(TS ts, PetscReal t, Vec U, Vec Udot, Vec F, void* ctx) {
  const PetscScalar *u, *udot;
  PetscScalar* f;
  VecGetArrayRead(U, &u);
  VecGetArrayRead(Udot, &udot);
  VecGetArray(F, &f);
  f[0] = udot[0] + u[0];
  f[1] = udot[1] - u[0];
  VecRestoreArrayRead(U, &u);
  VecRestoreArrayRead(Udot, &udot);
  VecRestoreArray(F, &f);
  return 0;
}

static PetscErrorCode IFunction_Nonconservative(TS ts, PetscReal t, Vec U, Vec Udot, Vec F, void* ctx) {
  const PetscScalar *u, *udot;
  PetscScalar* f;
  VecGetArrayRead(U, &u);
  VecGetArrayRead(Udot, &udot);
  VecGetArray(F, &f);
  f[0] = PetscExpScalar(u[0]) * udot[0] + PetscExpScalar(u[0]);
  f[1] = PetscExpScalar(u[1]) * udot[1] - PetscExpScalar(u[0]);
  VecRestoreArrayRead(U, &u);
  VecRestoreArrayRead(Udot, &udot);
  VecRestoreArray(F, &f);
  return 0;
}

static PetscErrorCode IFunction_TransientVar(TS ts, PetscReal t, Vec U, Vec Cdot, Vec F, void* ctx) {
  const PetscScalar *u, *cdot;
  PetscScalar* f;
  VecGetArrayRead(U, &u);
  VecGetArrayRead(Cdot, &cdot);
  VecGetArray(F, &f);
  f[0] = cdot[0] + PetscExpScalar(u[0]);
  f[1] = cdot[1] - PetscExpScalar(u[0]);
  VecRestoreArrayRead(U, &u);
  VecRestoreArrayRead(Cdot, &cdot);
  VecRestoreArray(F, &f);
  return 0;
}

static PetscErrorCode TransientVar(TS ts, Vec U, Vec C, void* ctx) {
  VecCopy(U, C);
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

// declare functions
static PetscErrorCode FormRHS(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec, ArrayPack& arp);
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
        arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y];
        arp.wtd(x, y) = 0.;
      } else {
        double rech_count = add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);
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
  PetscInt its;                // iterations for convergence
  SNESConvergedReason reason;  // Check convergence
  TS ts;
  Vec U;
  PetscScalar sum;

  // compute any starting values needed for arrays
  set_starting_values(params, arp);

  std::cout << "line 180" << std::endl;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetType(ts, TSBDF);
  TSSetDM(ts, user_context.da);
  VecCreateSeq(PETSC_COMM_SELF, 2, &U);
  VecSetValue(U, 0, 2., INSERT_VALUES);
  VecSetValue(U, 1, 1., INSERT_VALUES);
  std::cout << "line 187" << std::endl;
  VecLog(U);
  DMTSSetIFunction(user_context.da, IFunction_TransientVar, NULL);
  DMTSSetTransientVariable(user_context.da, TransientVar, NULL);
  std::cout << "line 191" << std::endl;
  TSSetMaxTime(ts, 1.);
  TSSetFromOptions(ts);
  TSSolve(ts, U);
  std::cout << "line 195" << std::endl;
  VecExp(U);

  VecView(U, PETSC_VIEWER_STDOUT_WORLD);
  VecSum(U, &sum);
  PetscPrintf(PETSC_COMM_WORLD, "Conservation error %g\n", (double)PetscRealPart(sum - 3.));
  VecDestroy(&U);
  TSDestroy(&ts);
  std::cout << "line 203" << std::endl;

  return 0;

  //
  //  // Get local array bounds
  //  const auto [xs, ys, xm, ym] = get_corners(user_context.da);
  //
  //  // values for storativity are reset each time; and wtd and recharge change from one timestep to the next, so set
  //  these
  //  // here
  //  for (auto j = ys; j < ys + ym; j++) {
  //    for (auto i = xs; i < xs + xm; i++) {
  //      dmdapack.S[i][j]        = arp.effective_storativity(i, j);
  //      dmdapack.h[i][j]        = arp.wtd(i, j);
  //      dmdapack.rech_vec[i][j] = arp.rech(i, j);
  //    }
  //  }
  //
  //  // Set local function evaluation routine
  //  DMDASNESSetFunctionLocal(
  //      user_context.da,
  //      INSERT_VALUES,
  //      (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal,
  //      &user_context);
  //
  //  // Evaluate initial guess
  //  FormInitialGuess(&user_context, user_context.da, user_context.x, arp);
  //
  //  // set the RHS
  //  FormRHS(&user_context, user_context.da, user_context.b, arp);
  //
  //  // Solve nonlinear system
  //  SNESSolve(user_context.snes, user_context.b, user_context.x);
  //  SNESGetIterationNumber(user_context.snes, &its);
  //  SNESGetConvergedReason(user_context.snes, &reason);
  //
  //  PetscPrintf(
  //      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason],
  //      its);
  //
  //  // recalculate the new storativity using the initial solve as an estimator for the final answer
  //  for (auto j = ys; j < ys + ym; j++) {
  //    for (auto i = xs; i < xs + xm; i++) {
  //      dmdapack.S[i][j] = updateEffectiveStorativity(
  //          arp.wtd(i, j), dmdapack.x[i][j] - arp.topo(i, j), arp.porosity(i, j), arp.effective_storativity(i, j));
  //    }
  //  }
  //
  //  // repeat the solve
  //  SNESSolve(user_context.snes, user_context.b, user_context.x);
  //  SNESGetIterationNumber(user_context.snes, &its);
  //  SNESGetConvergedReason(user_context.snes, &reason);
  //
  //  PetscPrintf(
  //      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason],
  //      its);
  //
  //  // copy the result back into the wtd array
  //  for (int j = ys; j < ys + ym; j++) {
  //    for (int i = xs; i < xs + xm; i++) {
  //      arp.wtd(i, j) = dmdapack.x[i][j] - arp.topo(i, j);
  //      if (arp.land_mask(i, j) == 0.f) {
  //        arp.total_loss_to_ocean +=
  //            arp.wtd(i, j) * arp.cell_area[j];  // could it be that because ocean cells are just set = x in the
  //            formula,
  //                                               // that loss to/gain from ocean is not properly recorded?
  //        arp.wtd(i, j) = 0.;
  //      }
  //    }
  //  }
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
