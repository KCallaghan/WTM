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
struct AppCtx {
  PetscReal timestep;
  PetscReal cellsize_NS;
  SNES snes       = nullptr;
  DM da           = nullptr;
  Vec x           = nullptr;  // Solution vector
  Vec b           = nullptr;  // RHS vector
  Vec S           = nullptr;
  Vec cellsize_EW = nullptr;
  Vec fdepth_vec  = nullptr;
  Vec ksat_vec    = nullptr;
  Vec mask        = nullptr;
  Vec porosity    = nullptr;
  Vec h           = nullptr;
  Vec topo_vec    = nullptr;

  ~AppCtx() {
    SNESDestroy(&snes);
    DMDestroy(&da);
    VecDestroy(&x);
    VecDestroy(&b);
    VecDestroy(&S);
    VecDestroy(&cellsize_EW);
    VecDestroy(&fdepth_vec);
    VecDestroy(&ksat_vec);
    VecDestroy(&mask);
    VecDestroy(&porosity);
    VecDestroy(&h);
    VecDestroy(&topo_vec);
  }

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &b);
    VecDuplicate(x, &S);
    VecDuplicate(x, &cellsize_EW);
    VecDuplicate(x, &fdepth_vec);
    VecDuplicate(x, &ksat_vec);
    VecDuplicate(x, &mask);
    VecDuplicate(x, &porosity);
    VecDuplicate(x, &h);
    VecDuplicate(x, &topo_vec);
  }
};

struct DMDA_Array_Pack {
  PetscScalar** x           = nullptr;
  PetscScalar** S           = nullptr;
  PetscScalar** cellsize_EW = nullptr;
  PetscScalar** fdepth_vec  = nullptr;
  PetscScalar** ksat_vec    = nullptr;
  PetscScalar** mask        = nullptr;
  PetscScalar** porosity    = nullptr;
  PetscScalar** h           = nullptr;
  PetscScalar** topo_vec    = nullptr;
  const AppCtx* context     = nullptr;

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    PETSC_CHECK(DMDAVecGetArray(user.da, user.x, &x));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.S, &S));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.cellsize_EW, &cellsize_EW));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.fdepth_vec, &fdepth_vec));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.ksat_vec, &ksat_vec));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.mask, &mask));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.porosity, &porosity));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.h, &h));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.topo_vec, &topo_vec));
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->x, &x));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->S, &S));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->cellsize_EW, &cellsize_EW));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->fdepth_vec, &fdepth_vec));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->ksat_vec, &ksat_vec));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->mask, &mask));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->porosity, &porosity));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->h, &h));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->topo_vec, &topo_vec));
    context = nullptr;
  }
};

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
        // Apply the first half of the recharge to the water-table depth grid (wtd)
        // use regular porosity for adding recharge since this checks
        // for underground space within add_recharge.
        arp.wtd(x, y) += add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);
      }
    }
  }
}

int update(Parameters& params, ArrayPack& arp) {
  AppCtx user;                /* user-defined work context */
  PetscInt its;               /* iterations for convergence */
  SNESConvergedReason reason; /* Check convergence */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize problem parameters
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.cellsize_NS = params.cellsize_n_s_metres;
  user.timestep    = params.deltat;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create nonlinear solver context
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  set_starting_values(params, arp);

  SNESCreate(PETSC_COMM_WORLD, &user.snes);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      DMDA_STENCIL_STAR,
      params.ncells_x,
      params.ncells_y,
      PETSC_DECIDE,
      PETSC_DECIDE,
      1,
      1,
      nullptr,
      nullptr,
      &user.da);
  DMSetFromOptions(user.da);
  DMSetUp(user.da);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DM; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.make_global_vectors();

  DMDA_Array_Pack dmdapack(user);

  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user.da);

  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.cellsize_EW[i][j] = arp.cellsize_e_w_metres[j];
      dmdapack.mask[i][j]        = arp.land_mask(i, j);
      dmdapack.fdepth_vec[i][j]  = arp.fdepth(i, j);
      dmdapack.ksat_vec[i][j]    = arp.ksat(i, j);
      dmdapack.S[i][j]           = arp.effective_storativity(i, j);
      dmdapack.porosity[i][j]    = arp.porosity(i, j);
      dmdapack.h[i][j]           = arp.wtd(i, j);
      dmdapack.topo_vec[i][j]    = arp.topo(i, j);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     User can override with:
     -snes_mf : matrix-free Newton-Krylov method with no preconditioning
                (unless user explicitly sets preconditioner)
     -snes_mf_operator : form preconditioning matrix as set by the user,
                         but use matrix-free approx for Jacobian-vector
                         products within Newton-Krylov method

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set local function evaluation routine
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMSetApplicationContext(user.da, &user);
  SNESSetDM(user.snes, user.da);

  DMDASNESSetFunctionLocal(
      user.da, INSERT_VALUES, (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal, &user);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESSetFromOptions(user.snes);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */
  FormInitialGuess(&user, user.da, user.x, arp);
  FormRHS(&user, user.da, user.b, arp);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESSolve(user.snes, user.b, user.x);
  SNESGetIterationNumber(user.snes, &its);
  SNESGetConvergedReason(user.snes, &reason);

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

  dmdapack.release();

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
static PetscErrorCode FormInitialGuess(AppCtx* user, DM da, Vec X, ArrayPack& arp) {
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
static PetscErrorCode FormRHS(AppCtx* user, DM da, Vec B, ArrayPack& arp) {
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
        b[i][j] = arp.topo(i, j) + arp.wtd(i, j);
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
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user) {
  DM da = user->da;
  PetscScalar **starting_storativity, **cellsize_ew, **my_mask, **my_fdepth, **my_ksat, **my_h, **my_porosity,
      **my_topo;
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
  PetscCall(DMDAVecGetArray(da, user->mask, &my_mask));
  PetscCall(DMDAVecGetArray(da, user->S, &starting_storativity));
  PetscCall(DMDAVecGetArray(da, user->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecGetArray(da, user->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecGetArray(da, user->ksat_vec, &my_ksat));
  PetscCall(DMDAVecGetArray(da, user->porosity, &my_porosity));
  PetscCall(DMDAVecGetArray(da, user->h, &my_h));
  PetscCall(DMDAVecGetArray(da, user->topo_vec, &my_topo));

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
            2. / (1. / depthIntegratedTransmissivity(x[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(x[i][j + 1], my_fdepth[i][j + 1], my_ksat[i][j + 1]));
        const PetscScalar e_W =
            2. / (1. / depthIntegratedTransmissivity(x[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(x[i][j - 1], my_fdepth[i][j - 1], my_ksat[i][j - 1]));
        const PetscScalar e_N =
            2. / (1. / depthIntegratedTransmissivity(x[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(x[i + 1][j], my_fdepth[i + 1][j], my_ksat[i + 1][j]));
        const PetscScalar e_S =
            2. / (1. / depthIntegratedTransmissivity(x[i][j], my_fdepth[i][j], my_ksat[i][j]) +
                  1. / depthIntegratedTransmissivity(x[i - 1][j], my_fdepth[i - 1][j], my_ksat[i - 1][j]));

        const PetscScalar uxx = -1. / (cellsize_ew[i][j] * cellsize_ew[i][j]) * (e_E * ux_E - e_W * ux_W);
        const PetscScalar uyy = -1. / (user->cellsize_NS * user->cellsize_NS) * (e_N * uy_N - e_S * uy_S);

        const PetscScalar my_S = updateEffectiveStorativity(
            my_h[i][j], x[i][j] - my_topo[i][j], my_porosity[i][j], starting_storativity[i][j]);

        f[i][j] = (uxx + uyy) * (user->timestep / my_S) + u;
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(da, user->mask, &my_mask));
  PetscCall(DMDAVecRestoreArray(da, user->S, &starting_storativity));
  PetscCall(DMDAVecRestoreArray(da, user->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecRestoreArray(da, user->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecRestoreArray(da, user->ksat_vec, &my_ksat));
  PetscCall(DMDAVecRestoreArray(da, user->porosity, &my_porosity));
  PetscCall(DMDAVecRestoreArray(da, user->h, &my_h));
  PetscCall(DMDAVecRestoreArray(da, user->topo_vec, &my_topo));

  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

}  // namespace FanDarcyGroundwater
