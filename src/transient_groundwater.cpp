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
  Vec r           = nullptr;  // Residual vector
  Vec b           = nullptr;  // RHS vector
  Vec T           = nullptr;
  Vec S           = nullptr;
  Vec cellsize_EW = nullptr;
  Vec mask        = nullptr;

  ~AppCtx() {
    SNESDestroy(&snes);
    DMDestroy(&da);
    VecDestroy(&x);
    VecDestroy(&r);
    VecDestroy(&b);
    VecDestroy(&T);
    VecDestroy(&S);
    VecDestroy(&cellsize_EW);
    VecDestroy(&mask);
  }

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &r);
    VecDuplicate(x, &b);
    VecDuplicate(x, &T);
    VecDuplicate(x, &S);
    VecDuplicate(x, &cellsize_EW);
    VecDuplicate(x, &mask);
  }
};

struct DMDA_Array_Pack {
  PetscScalar** x           = nullptr;
  PetscScalar** T           = nullptr;
  PetscScalar** S           = nullptr;
  PetscScalar** cellsize_EW = nullptr;
  PetscScalar** mask        = nullptr;
  const AppCtx* context     = nullptr;

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    PETSC_CHECK(DMDAVecGetArray(user.da, user.x, &x));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.T, &T));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.S, &S));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.cellsize_EW, &cellsize_EW));
    PETSC_CHECK(DMDAVecGetArray(user.da, user.mask, &mask));
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->x, &x));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->T, &T));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->S, &S));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->cellsize_EW, &cellsize_EW));
    PETSC_CHECK(DMDAVecRestoreArray(context->da, context->mask, &mask));
    context = nullptr;
  }
};

/*
   User-defined routines
*/
static PetscErrorCode FormRHS(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);
static PetscErrorCode FormJacobianLocal(DMDALocalInfo*, PetscScalar**, Mat, Mat, AppCtx*);

// User-defined routines
// static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);
// static PetscErrorCode FormJacobianLocal(DMDALocalInfo*, PetscScalar**, Mat, Mat, AppCtx*);
// static PetscErrorCode NonlinearGS(SNES, Vec, Vec, void*);

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
  // This assumes a ksat of 0.00005, which is an approximate global mean value, and an e-folding depth of 60.
  // e-folding depth of 60 corresponds to fdepth_a of 100, fdepth_b of 150, and slope of ~0.00443.
  // this should also be a reasonable e-folding depth for a range of other fdepth_a and _b parameters.
  // 1.5 is a standard value based on the shallow depths to which soil textures are known.
  constexpr double ocean_T = 0.00005 * (1.5 + 60.);

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
      if (arp.land_mask(x, y) == 0.f) {
        // in the ocean, we set several arrays to default values
        arp.transmissivity(x, y)        = ocean_T;
        arp.effective_storativity(x, y) = 1.;
      } else {
        // set the starting effective storativity
        if (arp.wtd(x, y) > 0) {
          arp.effective_storativity(x, y) = 1;
        } else {
          arp.effective_storativity(x, y) = arp.porosity(x, y);
        }
        // Apply the first half of the recharge to the water-table depth grid (wtd)
        // use regular porosity for adding recharge since this checks
        // for underground space within add_recharge.
        arp.wtd(x, y) += add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);

        arp.transmissivity(x, y) = depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y));
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

  // copy the transmissivity back into the T and storativity into the S vector
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.cellsize_EW[j][i] = arp.cellsize_e_w_metres[i];
      dmdapack.mask[j][i]        = arp.land_mask(j, i);
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
  // DMDASNESSetJacobianLocal(
  //   user.da, (PetscErrorCode(*)(DMDALocalInfo*, void*, Mat, Mat, void*))FormJacobianLocal, &user);

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

  for (int test = 0; test < 10; test++) {
    for (int y = 0; y < params.ncells_y; y++) {
      for (int x = 0; x < params.ncells_x; x++) {
        if (arp.land_mask(x, y) == 0.f) {
          // in the ocean, we set several arrays to default values
          arp.transmissivity(x, y) = 0.00005 * (1.5 + 60.);
        } else {
          arp.transmissivity(x, y) = depthIntegratedTransmissivity(dmdapack.x[x][y], arp.fdepth(x, y), arp.ksat(x, y));
          arp.effective_storativity(x, y) = updateEffectiveStorativity(
              arp.wtd(x, y), dmdapack.x[x][y], arp.porosity(x, y), arp.effective_storativity(x, y));
        }
      }
    }

    for (auto j = ys; j < ys + ym; j++) {
      for (auto i = xs; i < xs + xm; i++) {
        dmdapack.T[j][i] = arp.transmissivity(j, i);
        dmdapack.S[j][i] = arp.effective_storativity(j, i);
      }
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Solve nonlinear system
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    SNESSolve(user.snes, user.b, user.x);
    SNESGetIterationNumber(user.snes, &its);
    SNESGetConvergedReason(user.snes, &reason);

    PetscPrintf(
        PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  // copy the result back into the wtd array
  for (int j = ys; j < ys + ym; j++) {
    for (int i = xs; i < xs + xm; i++) {
      arp.wtd(j, i) = dmdapack.x[j][i] - arp.topo(j, i);
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
      if (arp.land_mask(j, i) == 0) {
        /* boundary conditions are all zero Dirichlet */
        x[j][i] = arp.topo(j, i) + 0.0;
      } else {
        x[j][i] = arp.topo(j, i) + arp.wtd(j, i);  // 0.0;
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
      if (arp.land_mask(j, i) == 0) {
        b[j][i] = arp.topo(j, i) + 0.0;
      } else {
        b[j][i] = arp.topo(j, i) + arp.wtd(j, i);  // 0.0;
      }
    }
  }
  DMDAVecRestoreArray(da, B, &b);
  return 0;
}

static inline PetscScalar eta(const AppCtx* ctx, int x, int y, int xplus, int yplus) {
  PetscScalar** my_T;

  PetscCall(DMDAVecGetArray(ctx->da, ctx->T, &my_T));

  const PetscScalar edge_T = 2. / (1. / my_T[x][y] + 1. / my_T[x + xplus][y + yplus]);

  PetscCall(DMDAVecRestoreArray(ctx->da, ctx->T, &my_T));

  return edge_T;
}
/* ------------------------------------------------------------------- */
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x).
 */
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user) {
  DM da = user->da;
  PetscScalar **my_S, **cellsize_ew, **my_mask;
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
  PetscCall(DMDAVecGetArray(da, user->S, &my_S));
  PetscCall(DMDAVecGetArray(da, user->cellsize_EW, &cellsize_ew));

  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      const PetscScalar u = x[j][i];
      if (my_mask[j][i] == 0) {
        f[j][i] = u;
      } else {
        const PetscScalar ux_E = (x[j][i + 1] - x[j][i]);
        const PetscScalar ux_W = (x[j][i] - x[j][i - 1]);
        const PetscScalar uy_N = (x[j + 1][i] - x[j][i]);
        const PetscScalar uy_S = (x[j][i] - x[j - 1][i]);
        const PetscScalar e_E  = eta(user, j, i, 0, 1);
        const PetscScalar e_W  = eta(user, j, i, 0, -1);
        const PetscScalar e_N  = eta(user, j, i, 1, 0);
        const PetscScalar e_S  = eta(user, j, i, -1, 0);
        const PetscScalar uxx  = -1. / (cellsize_ew[j][i] * cellsize_ew[j][i]) * (e_E * ux_E - e_W * ux_W);
        const PetscScalar uyy  = -1. / (user->cellsize_NS * user->cellsize_NS) * (e_N * uy_N - e_S * uy_S);

        /* For p=2, these terms decay to:
         uxx = (2.0*u - x[j][i-1] - x[j][i+1])*hydhx
         uyy = (2.0*u - x[j-1][i] - x[j+1][i])*hxdhy
        */
        f[j][i] = (uxx + uyy) * (user->timestep / my_S[j][i]) + u;
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(da, user->S, &my_S));
  PetscCall(DMDAVecRestoreArray(da, user->cellsize_EW, &cellsize_ew));
  PetscCall(DMDAVecRestoreArray(da, user->mask, &my_mask));

  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

/*
   FormJacobianLocal - Evaluates Jacobian matrix.
*/
static PetscErrorCode FormJacobianLocal(DMDALocalInfo* info, PetscScalar** x, Mat J, Mat B, AppCtx* user) {
  MatZeroEntries(B);
  PetscScalar **my_mask, **cellsize_ew;

  // Compute entries for the locally owned part of the Jacobian.
  //  - PETSc parallel matrix formats are partitioned by
  //    contiguous chunks of rows across the processors.
  //  - Each processor needs to insert only elements that it owns
  //    locally (but any non-local elements will be sent to the
  //    appropriate processor during matrix assembly).
  //  - Here, we set all entries for a particular row at once.
  PETSC_CHECK(DMDAVecGetArray(user->da, user->mask, &my_mask));
  PETSC_CHECK(DMDAVecGetArray(user->da, user->cellsize_EW, &cellsize_ew));

  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      MatStencil row{.k = 0, .j = j, .i = i, .c = 0};

      // boundary points
      if (my_mask[j][i] == 0) {
        std::array<PetscScalar, 1> local_v = {1.0};
        MatSetValuesStencil(B, 1, &row, 1, &row, local_v.data(), INSERT_VALUES);
      } else {
        // interior grid points
        // Jacobian from p=2
        std::array<MatStencil, 5> col;
        std::array<PetscScalar, 5> local_v;

        local_v[0] = 1.0 / (user->cellsize_NS * user->cellsize_NS);
        col[0].j   = j - 1;
        col[0].i   = i;
        local_v[1] = 1.0 / (cellsize_ew[j][i] * cellsize_ew[j][i]);
        col[1].j   = j;
        col[1].i   = i - 1;
        local_v[2] =
            2.0 * (1.0 / (user->cellsize_NS * user->cellsize_NS) + 1.0 / (cellsize_ew[j][i] * cellsize_ew[j][i]));

        col[2].j   = row.j;
        col[2].i   = row.i;
        local_v[3] = 1.0 / (cellsize_ew[j][i] * cellsize_ew[j][i]);
        col[3].j   = j;
        col[3].i   = i + 1;
        local_v[4] = 1.0 / (user->cellsize_NS * user->cellsize_NS);
        col[4].j   = j + 1;
        col[4].i   = i;
        MatSetValuesStencil(B, 1, &row, 5, col.data(), local_v.data(), INSERT_VALUES);
      }
    }
  }
  PETSC_CHECK(DMDAVecRestoreArray(user->da, user->mask, &my_mask));
  PETSC_CHECK(DMDAVecRestoreArray(user->da, user->cellsize_EW, &cellsize_ew));

  // Assemble matrix, using the 2-step process:
  //   MatAssemblyBegin(), MatAssemblyEnd().
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  if (J != B) {
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  }

  return 0;
}

}  // namespace FanDarcyGroundwater
