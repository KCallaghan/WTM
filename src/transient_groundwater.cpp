#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <omp.h>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;  // declares a row-major sparse matrix type of double

constexpr double seconds_in_a_year = 31536000.;

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal().
*/
typedef struct {
  PetscReal lambda, timestep; /* Bratu parameter */
  DM da;
  Vec T, S, b;
} AppCtx;
/*
   User-defined routines
*/
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
  SNES snes;                    /* nonlinear solver */
  Vec x, r, b, T, S;            /* solution, residual, rhs vectors */
  AppCtx user;                  /* user-defined work context */
  PetscInt its, xs, ys, xm, ym; /* iterations for convergence */
  SNESConvergedReason reason;   /* Check convergence */
  PetscBool alloc_star;         /* Only allocate for the STAR stencil  */
  SNESLineSearch linesearch;
  PetscScalar **my_T, **my_x, **my_S;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.lambda   = 0.0;
  user.timestep = params.deltat;
  alloc_star    = PETSC_FALSE;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "p-Bratu options", __FILE__);
  { PetscOptionsReal("-lambda", "Bratu parameter", "", user.lambda, &user.lambda, NULL); }
  PetscOptionsEnd();

  set_starting_values(params, arp);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      alloc_star ? DMDA_STENCIL_STAR : DMDA_STENCIL_BOX,
      params.ncells_x,
      params.ncells_y,
      PETSC_DECIDE,
      PETSC_DECIDE,
      1,
      1,
      NULL,
      NULL,
      &user.da);
  DMSetFromOptions(user.da);
  DMSetUp(user.da);
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DM; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(user.da, &x);
  VecDuplicate(x, &r);
  VecDuplicate(x, &b);
  user.b = b;
  PetscCall(VecDuplicate(x, &T));
  user.T = T;
  PetscCall(VecDuplicate(x, &S));
  user.S = S;

  PetscCall(DMDAVecGetArray(user.da, T, &my_T));
  PetscCall(DMDAVecGetArray(user.da, S, &my_S));

  PetscCall(DMDAGetCorners(user.da, &xs, &ys, NULL, &xm, &ym, NULL));

  // copy the transmissivity back into the T and storativity into the S vector
  for (int j = ys; j < ys + ym; j++)
    for (int i = xs; i < xs + xm; i++) {
      my_T[j][i] = arp.transmissivity(j, i);
      my_S[j][i] = arp.effective_storativity(j, i);
    }

  PetscCall(DMDAVecRestoreArray(user.da, T, &my_T));
  PetscCall(DMDAVecRestoreArray(user.da, S, &my_S));
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
  SNESSetDM(snes, user.da);

  DMDASNESSetFunctionLocal(
      user.da, INSERT_VALUES, (PetscErrorCode(*)(DMDALocalInfo*, void*, void*, void*))FormFunctionLocal, &user);
  DMDASNESSetJacobianLocal(
      user.da, (PetscErrorCode(*)(DMDALocalInfo*, void*, Mat, Mat, void*))FormJacobianLocal, &user);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESSetFromOptions(snes);
  SNESSetNGS(snes, NonlinearGS, &user);
  SNESGetLineSearch(snes, &linesearch);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */
  FormInitialGuess(&user, user.da, x, arp);
  FormRHS(&user, user.da, b, arp);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESSolve(snes, b, x);
  SNESGetIterationNumber(snes, &its);
  SNESGetConvergedReason(snes, &reason);
  PetscPrintf(
      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(DMDAVecGetArray(user.da, x, &my_x));
  PetscCall(DMDAGetCorners(user.da, &xs, &ys, NULL, &xm, &ym, NULL));

  // copy the result back into the wtd array
  for (int j = ys; j < ys + ym; j++)
    for (int i = xs; i < xs + xm; i++) {
      arp.wtd(j, i) = my_x[j][i] - arp.topo(j, i);
    }
  PetscCall(DMDAVecRestoreArray(user.da, x, &my_x));

  VecDestroy(&x);
  VecDestroy(&r);
  VecDestroy(&b);
  VecDestroy(&T);
  VecDestroy(&S);

  SNESDestroy(&snes);
  DMDestroy(&user.da);
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
  PetscInt i, j, Mx, My, xs, ys, xm, ym;
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
     Get a pointer to vector data.
       - For default PETSc vectors, VecGetArray() returns a pointer to
         the data array.  Otherwise, the routine is implementation dependent.
       - You MUST call VecRestoreArray() when you no longer need access to
         the array.
  */
  DMDAVecGetArray(da, X, &x);
  /*
     Get local grid boundaries (for 2-dimensional DA):
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)
  */
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  /*
     Compute initial guess over the locally owned part of the grid
  */
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      if (arp.land_mask(j, i) == 0) {
        /* boundary conditions are all zero Dirichlet */
        x[j][i] = 0.0;
      } else {
        x[j][i] = arp.wtd(j, i) + arp.topo(j, i);
      }
    }
  }
  /*
     Restore vector
  */
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
  PetscInt i, j, Mx, My, xs, ys, xm, ym;
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
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
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

/* p-Laplacian diffusivity */
// static inline PetscScalar eta(const AppCtx *ctx,int x,int y,int xplus,int yplus)
static inline PetscScalar eta(const AppCtx* ctx, int x, int y, PetscScalar ux, PetscScalar uy, int xplus, int yplus) {
  PetscInt xs, ys, xm, ym;
  DM da = ctx->da;
  PetscScalar** my_T;

  PetscCall(DMDAVecGetArray(da, ctx->T, &my_T));

  double edge_T = 2 / (1 / my_T[x][y] + 1 / my_T[x + xplus][y + yplus]);

  PetscCall(DMDAVecRestoreArray(da, ctx->T, &my_T));

  return edge_T;
}
static inline PetscScalar deta(const AppCtx* ctx, PetscReal x, PetscReal y, PetscScalar ux, PetscScalar uy) {
  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x).
 */
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user) {
  PetscReal hx, hy, dhx, dhy, sc;
  PetscInt i, j;
  DM da = user->da;

  PetscScalar **my_S, **my_b;
  hx  = 1.0 / (PetscReal)(info->mx - 1);
  hy  = 1.0 / (PetscReal)(info->my - 1);
  sc  = hx * hy * user->lambda;
  dhx = 1 / hx;
  dhy = 1 / hy;
  /*
     Compute function over the locally owned part of the grid
  */
  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      PetscReal xx = i * hx, yy = j * hy;
      if (i == 0 || j == 0 || i == info->mx - 1 || j == info->my - 1) {
        f[j][i] = x[j][i];
      } else {
        const PetscScalar u = x[j][i], ux_E = dhx * (x[j][i + 1] - x[j][i]),
                          uy_E = 0.25 * dhy * (x[j + 1][i] + x[j + 1][i + 1] - x[j - 1][i] - x[j - 1][i + 1]),
                          ux_W = dhx * (x[j][i] - x[j][i - 1]),
                          uy_W = 0.25 * dhy * (x[j + 1][i - 1] + x[j + 1][i] - x[j - 1][i - 1] - x[j - 1][i]),
                          ux_N = 0.25 * dhx * (x[j][i + 1] + x[j + 1][i + 1] - x[j][i - 1] - x[j + 1][i - 1]),
                          uy_N = dhy * (x[j + 1][i] - x[j][i]),
                          ux_S = 0.25 * dhx * (x[j - 1][i + 1] + x[j][i + 1] - x[j - 1][i - 1] - x[j][i - 1]),
                          uy_S = dhy * (x[j][i] - x[j - 1][i]), e_E = eta(user, i, j, ux_E, uy_E, 0, 1),
                          e_W = eta(user, i, j, ux_W, uy_W, 0, -1), e_N = eta(user, i, j, ux_N, uy_N, 1, 0),
                          e_S = eta(user, i, j, ux_S, uy_S, -1, 0), uxx = -hy * (e_E * ux_E - e_W * ux_W),
                          uyy = -hx * (e_N * uy_N - e_S * uy_S);

        PetscCall(DMDAVecGetArray(da, user->S, &my_S));
        PetscCall(DMDAVecGetArray(da, user->b, &my_b));

        f[j][i] = uxx + uyy + my_S[j][i] * ((x[j][i] - my_b[j][i]) / user->timestep);

        PetscCall(DMDAVecRestoreArray(da, user->S, &my_S));
        PetscCall(DMDAVecRestoreArray(da, user->b, &my_b));
      }
    }
  }
  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

/*
   FormJacobianLocal - Evaluates Jacobian matrix.
*/
static PetscErrorCode FormJacobianLocal(DMDALocalInfo* info, PetscScalar** x, Mat J, Mat B, AppCtx* user) {
  PetscInt i, j;
  MatStencil col[9], row;
  PetscScalar v[9];
  PetscReal hx, hy, hxdhy, hydhx, sc;
  MatZeroEntries(B);
  hx    = 1.0 / (PetscReal)(info->mx - 1);
  hy    = 1.0 / (PetscReal)(info->my - 1);
  sc    = hx * hy * user->lambda;
  hxdhy = hx / hy;
  hydhx = hy / hx;
  /*
     Compute entries for the locally owned part of the Jacobian.
      - PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Here, we set all entries for a particular row at once.
  */
  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      PetscReal xx = i * hx, yy = j * hy;
      row.j = j;
      row.i = i;
      /* boundary points */
      if (i == 0 || j == 0 || i == info->mx - 1 || j == info->my - 1) {
        v[0] = 1.0;
        MatSetValuesStencil(B, 1, &row, 1, &row, v, INSERT_VALUES);
      } else {
        /* interior grid points */
        const PetscScalar u = x[j][i];
        /* interior grid points */
        /* Jacobian from p=2 */
        v[0]     = -hxdhy;
        col[0].j = j - 1;
        col[0].i = i;
        v[1]     = -hydhx;
        col[1].j = j;
        col[1].i = i - 1;
        v[2]     = 2.0 * (hydhx + hxdhy) - sc * PetscExpScalar(u);
        col[2].j = row.j;
        col[2].i = row.i;
        v[3]     = -hydhx;
        col[3].j = j;
        col[3].i = i + 1;
        v[4]     = -hxdhy;
        col[4].j = j + 1;
        col[4].i = i;
        MatSetValuesStencil(B, 1, &row, 5, col, v, INSERT_VALUES);
      }
    }
  }
  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
  */
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
  PetscInt i, j, k, xs, ys, xm, ym, its, tot_its, sweeps, l, m;
  PetscReal hx, hy, hxdhy, hydhx, dhx, dhy, sc;
  PetscScalar **x, **b, bij, F, F0 = 0, J, y, u, eu;
  PetscReal atol, rtol, stol;
  DM da;
  AppCtx* user = (AppCtx*)ctx;
  Vec localX, localB;
  DMDALocalInfo info;
  SNESGetDM(snes, &da);
  DMDAGetLocalInfo(da, &info);
  hx    = 1.0 / (PetscReal)(info.mx - 1);
  hy    = 1.0 / (PetscReal)(info.my - 1);
  sc    = hx * hy * user->lambda;
  dhx   = 1 / hx;
  dhy   = 1 / hy;
  hxdhy = hx / hy;
  hydhx = hy / hx;

  tot_its = 0;
  SNESNGSGetSweeps(snes, &sweeps);
  SNESNGSGetTolerances(snes, &atol, &rtol, &stol, &its);
  DMGetLocalVector(da, &localX);
  if (B) {
    DMGetLocalVector(da, &localB);
  }
  if (B) {
    DMGlobalToLocalBegin(da, B, INSERT_VALUES, localB);
    DMGlobalToLocalEnd(da, B, INSERT_VALUES, localB);
  }
  if (B)
    DMDAVecGetArrayRead(da, localB, &b);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX);
  DMDAVecGetArray(da, localX, &x);
  for (l = 0; l < sweeps; l++) {
    /*
     Get local grid boundaries (for 2-dimensional DMDA):
     xs, ys   - starting grid indices (no ghost points)
     xm, ym   - widths of local grid (no ghost points)
     */
    DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
    for (m = 0; m < 2; m++) {
      for (j = ys; j < ys + ym; j++) {
        for (i = xs + (m + j) % 2; i < xs + xm; i += 2) {
          PetscReal xx = i * hx, yy = j * hy;
          if (B)
            bij = b[j][i];
          else
            bij = 0.;
          if (i == 0 || j == 0 || i == info.mx - 1 || j == info.my - 1) {
            /* boundary conditions are all zero Dirichlet */
            x[j][i] = 0.0 + bij;
          } else {
            const PetscScalar u_E = x[j][i + 1], u_W = x[j][i - 1], u_N = x[j + 1][i], u_S = x[j - 1][i];
            const PetscScalar uy_E = 0.25 * dhy * (x[j + 1][i] + x[j + 1][i + 1] - x[j - 1][i] - x[j - 1][i + 1]),
                              uy_W = 0.25 * dhy * (x[j + 1][i - 1] + x[j + 1][i] - x[j - 1][i - 1] - x[j - 1][i]),
                              ux_N = 0.25 * dhx * (x[j][i + 1] + x[j + 1][i + 1] - x[j][i - 1] - x[j + 1][i - 1]),
                              ux_S = 0.25 * dhx * (x[j - 1][i + 1] + x[j][i + 1] - x[j - 1][i - 1] - x[j][i - 1]);
            u                      = x[j][i];
            for (k = 0; k < its; k++) {
              const PetscScalar ux_E = dhx * (u_E - u), ux_W = dhx * (u - u_W), uy_N = dhy * (u_N - u),
                                uy_S = dhy * (u - u_S), e_E = eta(user, i, j, ux_E, uy_E, 0, 1),
                                e_W = eta(user, i, j, ux_W, uy_W, 0, -1), e_N = eta(user, i, j, ux_N, uy_N, 1, 0),
                                e_S = eta(user, i, j, ux_S, uy_S, -1, 0), de_E = deta(user, xx, yy, ux_E, uy_E),
                                de_W = deta(user, xx, yy, ux_W, uy_W), de_N = deta(user, xx, yy, ux_N, uy_N),
                                de_S = deta(user, xx, yy, ux_S, uy_S), newt_E = e_E + de_E * PetscSqr(ux_E),
                                newt_W = e_W + de_W * PetscSqr(ux_W), newt_N = e_N + de_N * PetscSqr(uy_N),
                                newt_S = e_S + de_S * PetscSqr(uy_S), uxx = -hy * (e_E * ux_E - e_W * ux_W),
                                uyy = -hx * (e_N * uy_N - e_S * uy_S);
              if (sc)
                eu = PetscExpScalar(u);
              else
                eu = 0;
              F = uxx + uyy - sc * eu - bij;
              if (k == 0)
                F0 = F;
              J = hxdhy * (newt_N + newt_S) + hydhx * (newt_E + newt_W) - sc * eu;
              y = F / J;
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
