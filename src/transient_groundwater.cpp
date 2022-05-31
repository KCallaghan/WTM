#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <omp.h>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

std::tuple<PetscInt, PetscInt, PetscInt, PetscInt> get_corners(const DM da) {
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(da, &xs, &ys, nullptr, &xm, &ym, nullptr);
  return {xs, ys, xm, ym};
}

// User-defined application context - contains data needed by the
// application-provided call-back routines, FormJacobianLocal() and
// FormFunctionLocal().
struct AppCtx {
  PetscReal lambda;  // Bratu parameter
  PetscReal timestep;
  PetscReal cellsize_NS;
  SNES snes;
  DM da;
  Vec x           = nullptr;  // Solution vector
  Vec r           = nullptr;  // Residual vector
  Vec b           = nullptr;  // RHS vector
  Vec T           = nullptr;
  Vec S           = nullptr;
  Vec H_T         = nullptr;
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
    VecDestroy(&H_T);
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
    VecDuplicate(x, &H_T);
    VecDuplicate(x, &cellsize_EW);
    VecDuplicate(x, &mask);
  }
};

struct DMDA_Array_Pack {
  PetscScalar** x           = nullptr;
  PetscScalar** T           = nullptr;
  PetscScalar** S           = nullptr;
  PetscScalar** H_T         = nullptr;
  PetscScalar** cellsize_EW = nullptr;
  PetscScalar** mask        = nullptr;
  const AppCtx* context     = nullptr;

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    DMDAVecGetArray(user.da, user.x, &x);
    DMDAVecGetArray(user.da, user.T, &T);
    DMDAVecGetArray(user.da, user.S, &S);
    DMDAVecGetArray(user.da, user.H_T, &H_T);
    DMDAVecGetArray(user.da, user.cellsize_EW, &cellsize_EW);
    DMDAVecGetArray(user.da, user.mask, &mask);
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    DMDAVecRestoreArray(context->da, context->x, &x);
    DMDAVecRestoreArray(context->da, context->T, &T);
    DMDAVecRestoreArray(context->da, context->S, &S);
    DMDAVecRestoreArray(context->da, context->H_T, &H_T);
    DMDAVecRestoreArray(context->da, context->cellsize_EW, &cellsize_EW);
    DMDAVecRestoreArray(context->da, context->mask, &mask);
    context = nullptr;
  }
};

// User-defined routines
static PetscErrorCode FormRHS(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);
static PetscErrorCode FormJacobianLocal(DMDALocalInfo*, PetscScalar**, Mat, Mat, AppCtx*);
static PetscErrorCode NonlinearGS(SNES, Vec, Vec, void*);

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
  AppCtx user;  // User-defined work context
  PetscInt iteration_count;
  SNESConvergedReason convergence_reason;

  // Initialize problem parameters
  user.lambda      = 0.0;
  user.timestep    = params.deltat;
  user.cellsize_NS = params.cellsize_n_s_metres;

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "p-Bratu options", __FILE__);
  PetscOptionsReal("-lambda", "Bratu parameter", "", user.lambda, &user.lambda, NULL);
  PetscOptionsEnd();

  set_starting_values(params, arp);

  // Create nonlinear solver and context
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &user.snes));

  // Create distributed array (DMDA) to manage parallel grid and vectors
  PetscCall(DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      DMDA_STENCIL_STAR,  // https://petsc.org/main/docs/manualpages/DMDA/DMDA_STENCIL_STAR/
      params.ncells_x,
      params.ncells_y,
      PETSC_DECIDE,
      PETSC_DECIDE,
      1,  // dof per node
      1,  // stencil width
      nullptr,
      nullptr,
      &user.da));
  PetscCall(DMSetFromOptions(user.da));
  PetscCall(DMSetUp(user.da));

  user.make_global_vectors();

  // Get DMDA arrays from the state arrays
  DMDA_Array_Pack dmdapack(user);

  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user.da);

  // copy the transmissivity back into the T and storativity into the S vector
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.T[j][i]           = arp.transmissivity(j, i);
      dmdapack.S[j][i]           = arp.effective_storativity(j, i);
      dmdapack.H_T[j][i]         = arp.wtd(j, i) + arp.topo(j, i);  // the current time's head
      dmdapack.cellsize_EW[j][i] = arp.cellsize_e_w_metres[i];
      dmdapack.mask[j][i]        = arp.land_mask(j, i);
    }
  }

  // User can override with:
  // -snes_mf : matrix-free Newton-Krylov method with no preconditioning
  //            (unless user explicitly sets preconditioner)
  // -snes_mf_operator : form preconditioning matrix as set by the user,
  //                     but use matrix-free approx for Jacobian-vector
  //                     products within Newton-Krylov method

  // Set local function evaluation routine
  DMSetApplicationContext(user.da, &user);
  SNESSetDM(user.snes, user.da);

  DMDASNESSetFunctionLocal(
      user.da, INSERT_VALUES, (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal, &user);
  DMDASNESSetJacobianLocal(
      user.da, (PetscErrorCode(*)(DMDALocalInfo*, void*, Mat, Mat, void*))FormJacobianLocal, &user);

  // Customize nonlinear solver; set runtime options
  SNESSetFromOptions(user.snes);
  SNESSetNGS(user.snes, NonlinearGS, &user);

  SNESLineSearch linesearch;
  SNESGetLineSearch(user.snes, &linesearch);

  // Evaluate initial guess
  // Note: The user should initialize the vector, x, with the initial guess
  // for the nonlinear solver prior to calling SNESSolve().  In particular,
  // to employ an initial guess of zero, the user should explicitly set
  // this vector to zero by calling VecSet().
  FormInitialGuess(&user, user.da, user.x, arp);
  FormRHS(&user, user.da, user.b, arp);

  // Solve nonlinear system
  SNESSolve(user.snes, user.b, user.x);
  SNESGetIterationNumber(user.snes, &iteration_count);
  SNESGetConvergedReason(user.snes, &convergence_reason);
  PetscPrintf(
      PETSC_COMM_WORLD,
      "%s Number of nonlinear iterations = %" PetscInt_FMT "\n",
      SNESConvergedReasons[convergence_reason],
      iteration_count);

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
static PetscErrorCode FormInitialGuess(AppCtx* /*user*/, DM da, Vec X, ArrayPack& arp) {
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

  // Get a pointer to vector data.
  //   - For default PETSc vectors, VecGetArray() returns a pointer to
  //     the data array.  Otherwise, the routine is implementation dependent.
  //   - You MUST call VecRestoreArray() when you no longer need access to
  //     the array.
  DMDAVecGetArray(da, X, &x);

  // Get local grid boundaries (for 2-dimensional DA):
  //   xs, ys   - starting grid indices (no ghost points)
  //   xm, ym   - widths of local grid (no ghost points)
  const auto [xs, ys, xm, ym] = get_corners(da);

  // Compute initial guess over the locally owned part of the grid
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      if (arp.land_mask(j, i) == 0) {
        // boundary conditions are all zero Dirichlet
        x[j][i] = 0.0;
      } else {
        x[j][i] = arp.wtd(j, i) + arp.topo(j, i);
      }
    }
  }

  // Restore vector
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
static PetscErrorCode FormRHS(AppCtx* /*user*/, DM da, Vec B, ArrayPack& arp) {
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
        b[j][i] = 0.0;
      } else {          // trying a case where the RHS is always set to 0
        b[j][i] = 0.0;  // arp.wtd(j, i) + arp.topo(j, i);
      }
    }
  }

  DMDAVecRestoreArray(da, B, &b);

  return 0;
}

// return the transmissivity at the cell edge, using the harmonic mean
static inline PetscScalar cell_edge_transmissivity(const AppCtx* ctx, int x, int y, int xplus, int yplus) {
  PetscScalar** my_T;

  PetscCall(DMDAVecGetArray(ctx->da, ctx->T, &my_T));

  const PetscScalar edge_T = 2 / (1 / my_T[x][y] + 1 / my_T[x + xplus][y + yplus]);

  PetscCall(DMDAVecRestoreArray(ctx->da, ctx->T, &my_T));

  return edge_T;
}

static inline PetscScalar
deta(const AppCtx* /*ctx*/, PetscReal /*x*/, PetscReal /*y*/, PetscScalar /*ux*/, PetscScalar /*uy*/) {
  return 0;
}

// FormFunctionLocal - Evaluates nonlinear function, F(x).
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user) {
  DM da = user->da;

  PetscScalar **my_S, **h_t, **cellsize_ew, **my_mask;
  const auto hx  = 1.0 / static_cast<PetscReal>(info->mx - 1);
  const auto hy  = 1.0 / static_cast<PetscReal>(info->my - 1);
  const auto dhx = 1 / hx;
  const auto dhy = 1 / hy;

  PetscCall(DMDAVecGetArray(da, user->mask, &my_mask));
  // Compute function over the locally owned part of the grid
  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      if (my_mask[j][i] == 0) {
        f[j][i] = x[j][i];
      } else {
        PetscCall(DMDAVecGetArray(da, user->S, &my_S));
        PetscCall(DMDAVecGetArray(da, user->H_T, &h_t));
        PetscCall(DMDAVecGetArray(da, user->cellsize_EW, &cellsize_ew));

        const PetscScalar ux_E = dhx * (x[j][i + 1] - x[j][i]);
        const PetscScalar ux_W = dhx * (x[j][i] - x[j][i - 1]);
        const PetscScalar uy_N = dhy * (x[j + 1][i] - x[j][i]);
        const PetscScalar uy_S = dhy * (x[j][i] - x[j - 1][i]);
        const PetscScalar e_E  = cell_edge_transmissivity(user, i, j, 0, 1);
        const PetscScalar e_W  = cell_edge_transmissivity(user, i, j, 0, -1);
        const PetscScalar e_N  = cell_edge_transmissivity(user, i, j, 1, 0);
        const PetscScalar e_S  = cell_edge_transmissivity(user, i, j, -1, 0);
        const PetscScalar uxx  = -1. / (cellsize_ew[j][i] * cellsize_ew[j][i]) * (e_E * ux_E - e_W * ux_W);
        const PetscScalar uyy  = -1. / (user->cellsize_NS * user->cellsize_NS) * (e_N * uy_N - e_S * uy_S);

        f[j][i] = uxx + uyy + my_S[j][i] * ((x[j][i] - h_t[j][i]) / user->timestep);

        PetscCall(DMDAVecRestoreArray(da, user->S, &my_S));
        PetscCall(DMDAVecRestoreArray(da, user->H_T, &h_t));
        PetscCall(DMDAVecRestoreArray(da, user->cellsize_EW, &cellsize_ew));
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(da, user->mask, &my_mask));

  PetscLogFlops(72.0 * info->xm * info->ym);

  return 0;
}

/*
   FormJacobianLocal - Evaluates Jacobian matrix.
*/
static PetscErrorCode FormJacobianLocal(DMDALocalInfo* info, PetscScalar** x, Mat J, Mat B, AppCtx* user) {
  MatZeroEntries(B);
  const auto hx    = 1.0 / static_cast<PetscReal>(info->mx - 1);
  const auto hy    = 1.0 / static_cast<PetscReal>(info->my - 1);
  const auto sc    = hx * hy * user->lambda;
  const auto hxdhy = hx / hy;
  const auto hydhx = hy / hx;
  PetscScalar local_v[9];
  PetscScalar** my_mask;
  DM da = user->da;

  // Compute entries for the locally owned part of the Jacobian.
  //  - PETSc parallel matrix formats are partitioned by
  //    contiguous chunks of rows across the processors.
  //  - Each processor needs to insert only elements that it owns
  //    locally (but any non-local elements will be sent to the
  //    appropriate processor during matrix assembly).
  //  - Here, we set all entries for a particular row at once.
  PetscCall(DMDAVecGetArray(da, user->mask, &my_mask));

  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      MatStencil row{.k = 0, .j = j, .i = i, .c = 0};

      // boundary points
      if (my_mask[j][i] == 0) {
        local_v[0] = 1.0;
        MatSetValuesStencil(B, 1, &row, 1, &row, local_v, INSERT_VALUES);
      } else {
        // interior grid points
        const PetscScalar u = x[j][i];
        // interior grid points
        // Jacobian from p=2
        MatStencil col[9];

        local_v[0] = -hxdhy;
        col[0].j   = j - 1;
        col[0].i   = i;
        local_v[1] = -hydhx;
        col[1].j   = j;
        col[1].i   = i - 1;
        local_v[2] = 2.0 * (hydhx + hxdhy) - sc * PetscExpScalar(u);
        col[2].j   = row.j;
        col[2].i   = row.i;
        local_v[3] = -hydhx;
        col[3].j   = j;
        col[3].i   = i + 1;
        local_v[4] = -hxdhy;
        col[4].j   = j + 1;
        col[4].i   = i;
        MatSetValuesStencil(B, 1, &row, 5, col, local_v, INSERT_VALUES);
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(da, user->mask, &my_mask));

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

/*
      Applies some sweeps on nonlinear Gauss-Seidel on each process
 */
PetscErrorCode NonlinearGS(SNES snes, Vec X, Vec B, void* ctx) {
  PetscScalar **x, **b, bij, F0 = 0;
  DM da;
  AppCtx* const user = static_cast<AppCtx*>(ctx);
  Vec localX, localB;
  DMDALocalInfo info;
  SNESGetDM(snes, &da);
  DMDAGetLocalInfo(da, &info);
  const auto hx    = 1.0 / static_cast<PetscReal>(info.mx - 1);
  const auto hy    = 1.0 / static_cast<PetscReal>(info.my - 1);
  const auto sc    = hx * hy * user->lambda;
  const auto dhx   = 1 / hx;
  const auto dhy   = 1 / hy;
  const auto hxdhy = hx / hy;
  const auto hydhx = hy / hx;

  PetscInt tot_its = 0;

  PetscInt sweeps;
  SNESNGSGetSweeps(snes, &sweeps);

  PetscInt iterations;
  PetscReal atol, rtol, stol;
  SNESNGSGetTolerances(snes, &atol, &rtol, &stol, &iterations);

  DMGetLocalVector(da, &localX);
  if (B) {
    DMGetLocalVector(da, &localB);
  }
  if (B) {
    DMGlobalToLocalBegin(da, B, INSERT_VALUES, localB);
    DMGlobalToLocalEnd(da, B, INSERT_VALUES, localB);
  }
  if (B) {
    DMDAVecGetArrayRead(da, localB, &b);
  }
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
  DMDAVecGetArray(da, localX, &x);
  for (auto l = 0; l < sweeps; l++) {
    // Get local grid boundaries (for 2-dimensional DMDA):
    // xs, ys   - starting grid indices (no ghost points)
    // xm, ym   - widths of local grid (no ghost points)
    const auto [xs, ys, xm, ym] = get_corners(da);
    for (auto m = 0; m < 2; m++) {
      for (auto j = ys; j < ys + ym; j++) {
        for (auto i = xs + (m + j) % 2; i < xs + xm; i += 2) {
          PetscReal xx = i * hx, yy = j * hy;
          if (B) {
            bij = b[j][i];
          } else {
            bij = 0.;
          }
          if (i == 0 || j == 0 || i == info.mx - 1 || j == info.my - 1) {
            /* boundary conditions are all zero Dirichlet */
            x[j][i] = 0.0 + bij;
          } else {
            const PetscScalar u_E  = x[j][i + 1];
            const PetscScalar u_W  = x[j][i - 1];
            const PetscScalar u_N  = x[j + 1][i];
            const PetscScalar u_S  = x[j - 1][i];
            const PetscScalar uy_E = 0.25 * dhy * (x[j + 1][i] + x[j + 1][i + 1] - x[j - 1][i] - x[j - 1][i + 1]);
            const PetscScalar uy_W = 0.25 * dhy * (x[j + 1][i - 1] + x[j + 1][i] - x[j - 1][i - 1] - x[j - 1][i]);
            const PetscScalar ux_N = 0.25 * dhx * (x[j][i + 1] + x[j + 1][i + 1] - x[j][i - 1] - x[j + 1][i - 1]);
            const PetscScalar ux_S = 0.25 * dhx * (x[j - 1][i + 1] + x[j][i + 1] - x[j - 1][i - 1] - x[j][i - 1]);
            PetscScalar u          = x[j][i];
            for (auto k = 0; k < iterations; k++) {
              const PetscScalar ux_E   = dhx * (u_E - u);
              const PetscScalar ux_W   = dhx * (u - u_W);
              const PetscScalar uy_N   = dhy * (u_N - u);
              const PetscScalar uy_S   = dhy * (u - u_S);
              const PetscScalar e_E    = cell_edge_transmissivity(user, i, j, 0, 1);
              const PetscScalar e_W    = cell_edge_transmissivity(user, i, j, 0, -1);
              const PetscScalar e_N    = cell_edge_transmissivity(user, i, j, 1, 0);
              const PetscScalar e_S    = cell_edge_transmissivity(user, i, j, -1, 0);
              const PetscScalar de_E   = deta(user, xx, yy, ux_E, uy_E);
              const PetscScalar de_W   = deta(user, xx, yy, ux_W, uy_W);
              const PetscScalar de_N   = deta(user, xx, yy, ux_N, uy_N);
              const PetscScalar de_S   = deta(user, xx, yy, ux_S, uy_S);
              const PetscScalar newt_E = e_E + de_E * PetscSqr(ux_E);
              const PetscScalar newt_W = e_W + de_W * PetscSqr(ux_W);
              const PetscScalar newt_N = e_N + de_N * PetscSqr(uy_N);
              const PetscScalar newt_S = e_S + de_S * PetscSqr(uy_S);
              const PetscScalar uxx    = -hy * (e_E * ux_E - e_W * ux_W);
              const PetscScalar uyy    = -hx * (e_N * uy_N - e_S * uy_S);

              const PetscScalar eu = sc ? PetscExpScalar(u) : 0;

              const auto F = uxx + uyy - sc * eu - bij;
              if (k == 0) {
                F0 = F;
              }
              const auto J = hxdhy * (newt_N + newt_S) + hydhx * (newt_E + newt_W) - sc * eu;
              const auto y = F / J;
              u -= y;
              tot_its++;
              if (atol > PetscAbsReal(PetscRealPart(F)) ||
                  rtol * PetscAbsReal(PetscRealPart(F0)) > PetscAbsReal(PetscRealPart(F)) ||
                  stol * PetscAbsReal(PetscRealPart(u)) > PetscAbsReal(PetscRealPart(y))) {
                break;
              }
            }
            x[j][i] = u;
          }
        }
      }
    }
    /*
x     Restore vector
     */
  }
  DMDAVecRestoreArray(da, localX, &x);
  DMLocalToGlobalBegin(da, localX, INSERT_VALUES, X);
  DMLocalToGlobalEnd(da, localX, INSERT_VALUES, X);
  PetscLogFlops(tot_its * (118.0));
  DMRestoreLocalVector(da, &localX);
  if (B) {
    DMDAVecRestoreArrayRead(da, localB, &b);
    DMRestoreLocalVector(da, &localB);
  }

  return 0;
}

}  // namespace FanDarcyGroundwater
