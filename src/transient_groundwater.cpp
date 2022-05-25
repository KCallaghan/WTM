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

typedef enum { JAC_BRATU, JAC_PICARD, JAC_STAR, JAC_NEWTON } JacType;
static const char* const JacTypes[] = {"BRATU", "PICARD", "STAR", "NEWTON", "JacType", "JAC_", 0};
/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal().
*/
typedef struct {
  PetscReal lambda;  /* Bratu parameter */
  PetscReal p;       /* Exponent in p-Laplacian */
  PetscReal epsilon; /* Regularization */
  PetscReal source;  /* Source term */
  JacType jtype;     /* What type of Jacobian to assemble */
  PetscBool picard;
  PetscInt blocks[2];
  PetscReal kappa;
  PetscInt initial; /* initial conditions type */
  DM da;
  Vec T;
} AppCtx;
/*
   User-defined routines
*/
static PetscErrorCode FormRHS(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormInitialGuess(AppCtx*, DM, Vec, ArrayPack& arp);
static PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscScalar**, PetscScalar**, AppCtx*);
static PetscErrorCode FormJacobianLocal(DMDALocalInfo*, PetscScalar**, Mat, Mat, AppCtx*);
static PetscErrorCode NonlinearGS(SNES, Vec, Vec, void*);
typedef struct _n_PreCheck* PreCheck;
struct _n_PreCheck {
  MPI_Comm comm;
  PetscReal angle;
  Vec Ylast;
  PetscViewer monitor;
};
PetscErrorCode PreCheckCreate(MPI_Comm, PreCheck*);
PetscErrorCode PreCheckDestroy(PreCheck*);
PetscErrorCode PreCheckFunction(SNESLineSearch, Vec, Vec, PetscBool*, void*);
PetscErrorCode PreCheckSetFromOptions(PreCheck);

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
        arp.original_wtd(x, y)          = 0.;
        arp.effective_storativity(x, y) = 1.;
      } else {
        arp.original_wtd(x, y) = arp.wtd(x, y);
        // Apply the first half of the recharge to the water-table depth grid (wtd)
        // use regular porosity for adding recharge since this checks
        // for underground space within add_recharge.
        arp.wtd(x, y) += add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);

        arp.transmissivity(x, y) = depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y));
        // also set the starting effective storativity
        if (arp.original_wtd(x, y) > 0) {
          arp.effective_storativity(x, y) = 1;
        } else {
          arp.effective_storativity(x, y) = arp.porosity(x, y);
        }
      }
      // std::cout<<"x "<<x<<" y "<<y<<" transmissivity "<<arp.transmissivity(x,y)<<std::endl;
    }
  }
}

int update(Parameters& params, ArrayPack& arp) {
  SNES snes;                    /* nonlinear solver */
  Vec x, r, b, T;               /* solution, residual, rhs vectors */
  AppCtx user;                  /* user-defined work context */
  PetscInt its, xs, ys, xm, ym; /* iterations for convergence */
  SNESConvergedReason reason;   /* Check convergence */
  PetscBool alloc_star;         /* Only allocate for the STAR stencil  */
  PetscBool write_output;
  char filename[PETSC_MAX_PATH_LEN] = "ex15.vts";
  PetscReal bratu_lambda_max = 6.81, bratu_lambda_min = 0.;
  PreCheck precheck = NULL; /* precheck context for version in this file */
  PetscInt use_precheck;    /* 0=none, 1=version in this file, 2=SNES-provided version */
  PetscReal precheck_angle; /* When manually setting the SNES-provided precheck function */
  SNESLineSearch linesearch;
  PetscScalar **my_T, **my_x;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Initialize problem parameters
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.lambda    = 0.0;
  user.p         = 2.0;
  user.epsilon   = 1e-5;
  user.source    = 0.1;
  user.jtype     = JAC_NEWTON;
  user.initial   = -1;
  user.blocks[0] = 1;
  user.blocks[1] = 1;
  user.kappa     = 1e-3;
  alloc_star     = PETSC_FALSE;
  use_precheck   = 0;
  precheck_angle = 10.;
  user.picard    = PETSC_FALSE;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "p-Bratu options", __FILE__);
  {
    PetscInt two = 2;
    PetscOptionsReal("-lambda", "Bratu parameter", "", user.lambda, &user.lambda, NULL);
    PetscOptionsReal("-p", "Exponent `p' in p-Laplacian", "", user.p, &user.p, NULL);
    PetscOptionsReal("-epsilon", "Strain-regularization in p-Laplacian", "", user.epsilon, &user.epsilon, NULL);
    PetscOptionsReal("-source", "Constant source term", "", user.source, &user.source, NULL);
    PetscOptionsEnum(
        "-jtype",
        "Jacobian approximation to assemble",
        "",
        JacTypes,
        (PetscEnum)user.jtype,
        (PetscEnum*)&user.jtype,
        NULL);
    PetscOptionsName("-picard", "Solve with defect-correction Picard iteration", "", &user.picard);
    if (user.picard) {
      user.jtype = JAC_PICARD;
      /* the Picard linearization only requires a 5 point stencil, while the Newton linearization requires a 9 point
       * stencil */
      /* hence allocating the 5 point stencil gives the same convergence as the 9 point stencil since the extra stencil
       * points are not used */
      PetscOptionsBool("-alloc_star", "Allocate for STAR stencil (5-point)", "", alloc_star, &alloc_star, NULL);
    }
    PetscOptionsInt(
        "-precheck",
        "Use a pre-check correction intended for use with Picard iteration 1=this version, 2=library",
        "",
        use_precheck,
        &use_precheck,
        NULL);
    PetscOptionsInt(
        "-initial",
        "Initial conditions type (-1: default, 0: zero-valued, 1: peaked guess)",
        "",
        user.initial,
        &user.initial,
        NULL);
    if (use_precheck == 2) { /* Using library version, get the angle */
      PetscOptionsReal(
          "-precheck_angle",
          "Angle in degrees between successive search directions necessary to activate step correction",
          "",
          precheck_angle,
          &precheck_angle,
          NULL);
    }
    PetscOptionsIntArray(
        "-blocks", "number of coefficient interfaces in x and y direction", "", user.blocks, &two, NULL);
    PetscOptionsReal("-kappa", "diffusivity in odd regions", "", user.kappa, &user.kappa, NULL);
    PetscOptionsString("-o", "Output solution in vts format", "", filename, filename, sizeof(filename), &write_output);
  }
  PetscOptionsEnd();
  if (user.lambda > bratu_lambda_max || user.lambda < bratu_lambda_min) {
    PetscPrintf(PETSC_COMM_WORLD, "WARNING: lambda %g out of range for p=2\n", (double)user.lambda);
  }

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
      4,
      4,
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
  PetscCall(VecDuplicate(x, &T));
  user.T = T;

  PetscCall(DMDAVecGetArray(user.da, T, &my_T));
  PetscCall(DMDAGetCorners(user.da, &xs, &ys, NULL, &xm, &ym, NULL));

  // copy the transmissivity back into the T vector
  for (int j = ys; j < ys + ym; j++)
    for (int i = xs; i < xs + xm; i++) {
      my_T[j][i] = 1.;  // arp.transmissivity(j, i);
      std::cout << "j " << j << " i " << i << " transmissivity " << arp.transmissivity(j, i) << std::endl;
    }
  PetscCall(DMDAVecRestoreArray(user.da, T, &my_T));
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
  /* Set up the precheck context if requested */
  if (use_precheck == 1) { /* Use the precheck routines in this file */
    PreCheckCreate(PETSC_COMM_WORLD, &precheck);
    PreCheckSetFromOptions(precheck);
    SNESLineSearchSetPreCheck(linesearch, PreCheckFunction, precheck);
  } else if (use_precheck == 2) { /* Use the version provided by the library */
    SNESLineSearchSetPreCheck(linesearch, SNESLineSearchPreCheckPicard, &precheck_angle);
  }
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
  if (write_output) {
    PetscViewer viewer;
    PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
  }
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
  SNESDestroy(&snes);
  DMDestroy(&user.da);
  PreCheckDestroy(&precheck);
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
  PetscReal temp1, temp, hx, hy;
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
  hx    = 1.0 / (PetscReal)(Mx - 1);
  hy    = 1.0 / (PetscReal)(My - 1);
  temp1 = user->lambda / (user->lambda + 1.);
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
    temp = (PetscReal)(PetscMin(j, My - j - 1)) * hy;
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
  PetscReal hx, hy;
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
  hx = 1.0 / (PetscReal)(Mx - 1);
  hy = 1.0 / (PetscReal)(My - 1);
  DMDAVecGetArray(da, B, &b);
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      if (arp.land_mask(j, i) == 0) {
        b[j][i] = 0.0;
      } else {
        b[j][i] = arp.wtd(j, i) + arp.topo(j, i);  // hx*hy*user->source;
      }
    }
  }
  DMDAVecRestoreArray(da, B, &b);
  return 0;
}
static inline PetscReal kappa(const AppCtx* ctx, PetscReal x, PetscReal y) {
  return (((PetscInt)(x * ctx->blocks[0])) + ((PetscInt)(y * ctx->blocks[1]))) % 2 ? ctx->kappa : 1.0;
}
/* p-Laplacian diffusivity */
// static inline PetscScalar eta(const AppCtx *ctx,int x,int y,int xplus,int yplus)
static inline PetscScalar eta(const AppCtx* ctx, int x, int y, PetscScalar ux, PetscScalar uy) {
  PetscInt xs, ys, xm, ym;
  DM da = ctx->da;
  PetscScalar** my_T;
  int xplus = 0;
  int yplus = 0;
  PetscCall(DMDAVecGetArray(da, ctx->T, &my_T));
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  double edge_T = 2 / (1 / my_T[x][y] + 1 / my_T[x + xplus][y + yplus]);
  // std::cout<<"my T is "<<my_T[x][y];
  PetscCall(DMDAVecRestoreArray(da, ctx->T, &my_T));

  return edge_T;
  // return kappa(ctx,x,y) * PetscPowScalar(PetscSqr(ctx->epsilon)+0.5*(ux*ux + uy*uy),0.5*(ctx->p-2.));
}
static inline PetscScalar deta(const AppCtx* ctx, PetscReal x, PetscReal y, PetscScalar ux, PetscScalar uy) {
  return (ctx->p == 2) ? 0
                       : kappa(ctx, x, y) *
                             PetscPowScalar(PetscSqr(ctx->epsilon) + 0.5 * (ux * ux + uy * uy), 0.5 * (ctx->p - 4)) *
                             0.5 * (ctx->p - 2.);
}
/* ------------------------------------------------------------------- */
/*
   FormFunctionLocal - Evaluates nonlinear function, F(x).
 */
static PetscErrorCode FormFunctionLocal(DMDALocalInfo* info, PetscScalar** x, PetscScalar** f, AppCtx* user) {
  PetscReal hx, hy, dhx, dhy, sc;
  PetscInt i, j;
  PetscScalar eu;
  hx  = 1.0 / (PetscReal)(info->mx - 1);
  hy  = 1.0 / (PetscReal)(info->my - 1);
  sc  = hx * hy * user->lambda;
  dhx = 1 / hx;
  dhy = 1 / hy;
  /*
     Compute function over the locally owned part of the grid
  */
  // std::cout<<"mx "<<info->mx<<" my "<<info->my<<std::endl;
  for (j = info->ys; j < info->ys + info->ym; j++) {
    for (i = info->xs; i < info->xs + info->xm; i++) {
      //  std::cout<<"j "<<j<<" i "<<i<<std::endl;
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
                          uy_S = dhy * (x[j][i] - x[j - 1][i]), e_E = eta(user, i, j, ux_E, uy_E),
                          e_W = eta(user, i, j, ux_W, uy_W), e_N = eta(user, i, j, ux_N, uy_N),
                          e_S = eta(user, i, j, ux_S, uy_S), uxx = -hy * (e_E * ux_E - e_W * ux_W),
                          uyy = -hx * (e_N * uy_N - e_S * uy_S);
        std::cout << "i " << i << " j " << j << "e_E " << e_E << " ux_E " << ux_E << std::endl;
        if (sc)
          eu = PetscExpScalar(u);
        else
          eu = 0.;
        /* For p=2, these terms decay to:
         uxx = (2.0*u - x[j][i-1] - x[j][i+1])*hydhx
         uyy = (2.0*u - x[j-1][i] - x[j+1][i])*hxdhy
        */
        f[j][i] = uxx + uyy - sc * eu;
        std::cout << "f " << f[j][i] << std::endl;
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
  PetscReal hx, hy, hxdhy, hydhx, dhx, dhy, sc;
  MatZeroEntries(B);
  hx    = 1.0 / (PetscReal)(info->mx - 1);
  hy    = 1.0 / (PetscReal)(info->my - 1);
  sc    = hx * hy * user->lambda;
  dhx   = 1 / hx;
  dhy   = 1 / hy;
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
        const PetscScalar ux_E = dhx * (x[j][i + 1] - x[j][i]),
                          uy_E = 0.25 * dhy * (x[j + 1][i] + x[j + 1][i + 1] - x[j - 1][i] - x[j - 1][i + 1]),
                          ux_W = dhx * (x[j][i] - x[j][i - 1]),
                          uy_W = 0.25 * dhy * (x[j + 1][i - 1] + x[j + 1][i] - x[j - 1][i - 1] - x[j - 1][i]),
                          ux_N = 0.25 * dhx * (x[j][i + 1] + x[j + 1][i + 1] - x[j][i - 1] - x[j + 1][i - 1]),
                          uy_N = dhy * (x[j + 1][i] - x[j][i]),
                          ux_S = 0.25 * dhx * (x[j - 1][i + 1] + x[j][i + 1] - x[j - 1][i - 1] - x[j][i - 1]),
                          uy_S = dhy * (x[j][i] - x[j - 1][i]), u = x[j][i], e_E = eta(user, i, j, ux_E, uy_E),
                          e_W = eta(user, i, j, ux_W, uy_W), e_N = eta(user, i, j, ux_N, uy_N),
                          e_S = eta(user, i, j, ux_S, uy_S), de_E = deta(user, xx, yy, ux_E, uy_E),
                          de_W = deta(user, xx, yy, ux_W, uy_W), de_N = deta(user, xx, yy, ux_N, uy_N),
                          de_S = deta(user, xx, yy, ux_S, uy_S), skew_E = de_E * ux_E * uy_E,
                          skew_W = de_W * ux_W * uy_W, skew_N = de_N * ux_N * uy_N, skew_S = de_S * ux_S * uy_S,
                          cross_EW = 0.25 * (skew_E - skew_W), cross_NS = 0.25 * (skew_N - skew_S),
                          newt_E = e_E + de_E * PetscSqr(ux_E), newt_W = e_W + de_W * PetscSqr(ux_W),
                          newt_N = e_N + de_N * PetscSqr(uy_N), newt_S = e_S + de_S * PetscSqr(uy_S);
        /* interior grid points */
        switch (user->jtype) {
          case JAC_BRATU:
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
            break;
          case JAC_PICARD:
            /* Jacobian arising from Picard linearization */
            v[0]     = -hxdhy * e_S;
            col[0].j = j - 1;
            col[0].i = i;
            v[1]     = -hydhx * e_W;
            col[1].j = j;
            col[1].i = i - 1;
            v[2]     = (e_W + e_E) * hydhx + (e_S + e_N) * hxdhy;
            col[2].j = row.j;
            col[2].i = row.i;
            v[3]     = -hydhx * e_E;
            col[3].j = j;
            col[3].i = i + 1;
            v[4]     = -hxdhy * e_N;
            col[4].j = j + 1;
            col[4].i = i;
            MatSetValuesStencil(B, 1, &row, 5, col, v, INSERT_VALUES);
            break;
          case JAC_STAR:
            /* Full Jacobian, but only a star stencil */
            col[0].j = j - 1;
            col[0].i = i;
            col[1].j = j;
            col[1].i = i - 1;
            col[2].j = j;
            col[2].i = i;
            col[3].j = j;
            col[3].i = i + 1;
            col[4].j = j + 1;
            col[4].i = i;
            v[0]     = -hxdhy * newt_S + cross_EW;
            v[1]     = -hydhx * newt_W + cross_NS;
            v[2]     = hxdhy * (newt_N + newt_S) + hydhx * (newt_E + newt_W) - sc * PetscExpScalar(u);
            v[3]     = -hydhx * newt_E - cross_NS;
            v[4]     = -hxdhy * newt_N - cross_EW;
            MatSetValuesStencil(B, 1, &row, 5, col, v, INSERT_VALUES);
            break;
          case JAC_NEWTON:
            /* The Jacobian is
             -div [ eta (grad u) + deta (grad u0 . grad u) grad u0 ] - (eE u0) u
            */
            col[0].j = j - 1;
            col[0].i = i - 1;
            col[1].j = j - 1;
            col[1].i = i;
            col[2].j = j - 1;
            col[2].i = i + 1;
            col[3].j = j;
            col[3].i = i - 1;
            col[4].j = j;
            col[4].i = i;
            col[5].j = j;
            col[5].i = i + 1;
            col[6].j = j + 1;
            col[6].i = i - 1;
            col[7].j = j + 1;
            col[7].i = i;
            col[8].j = j + 1;
            col[8].i = i + 1;
            v[0]     = -0.25 * (skew_S + skew_W);
            v[1]     = -hxdhy * newt_S + cross_EW;
            v[2]     = 0.25 * (skew_S + skew_E);
            v[3]     = -hydhx * newt_W + cross_NS;
            v[4]     = hxdhy * (newt_N + newt_S) + hydhx * (newt_E + newt_W) - sc * PetscExpScalar(u);
            v[5]     = -hydhx * newt_E - cross_NS;
            v[6]     = 0.25 * (skew_N + skew_W);
            v[7]     = -hxdhy * newt_N - cross_EW;
            v[8]     = -0.25 * (skew_N + skew_E);
            MatSetValuesStencil(B, 1, &row, 9, col, v, INSERT_VALUES);
            break;
          default:
            SETERRQ(
                PetscObjectComm((PetscObject)info->da), PETSC_ERR_SUP, "Jacobian type %d not implemented", user->jtype);
        }
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
  /*
     Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.
  */
  if (user->jtype == JAC_NEWTON) {
    PetscLogFlops(info->xm * info->ym * (131.0));
  }
  return 0;
}
/***********************************************************
 * PreCheck implementation
 ***********************************************************/
PetscErrorCode PreCheckSetFromOptions(PreCheck precheck) {
  PetscBool flg;
  PetscOptionsBegin(precheck->comm, NULL, "PreCheck Options", "none");
  PetscOptionsReal(
      "-precheck_angle",
      "Angle in degrees between successive search directions necessary to activate step correction",
      "",
      precheck->angle,
      &precheck->angle,
      NULL);
  flg = PETSC_FALSE;
  PetscOptionsBool("-precheck_monitor", "Monitor choices made by precheck routine", "", flg, &flg, NULL);
  if (flg) {
    PetscViewerASCIIOpen(precheck->comm, "stdout", &precheck->monitor);
  }
  PetscOptionsEnd();
  return 0;
}
/*
  Compare the direction of the current and previous step, modify the current step accordingly
*/
PetscErrorCode PreCheckFunction(SNESLineSearch linesearch, Vec X, Vec Y, PetscBool* changed, void* ctx) {
  PreCheck precheck;
  Vec Ylast;
  PetscScalar dot;
  PetscInt iter;
  PetscReal ynorm, ylastnorm, theta, angle_radians;
  SNES snes;
  SNESLineSearchGetSNES(linesearch, &snes);
  precheck = (PreCheck)ctx;
  if (!precheck->Ylast)
    VecDuplicate(Y, &precheck->Ylast);
  Ylast = precheck->Ylast;
  SNESGetIterationNumber(snes, &iter);
  if (iter < 1) {
    VecCopy(Y, Ylast);
    *changed = PETSC_FALSE;
    return 0;
  }
  VecDot(Y, Ylast, &dot);
  VecNorm(Y, NORM_2, &ynorm);
  VecNorm(Ylast, NORM_2, &ylastnorm);
  /* Compute the angle between the vectors Y and Ylast, clip to keep inside the domain of acos() */
  theta         = PetscAcosReal((PetscReal)PetscClipInterval(PetscAbsScalar(dot) / (ynorm * ylastnorm), -1.0, 1.0));
  angle_radians = precheck->angle * PETSC_PI / 180.;
  if (PetscAbsReal(theta) < angle_radians || PetscAbsReal(theta - PETSC_PI) < angle_radians) {
    /* Modify the step Y */
    PetscReal alpha, ydiffnorm;
    VecAXPY(Ylast, -1.0, Y);
    VecNorm(Ylast, NORM_2, &ydiffnorm);
    alpha = ylastnorm / ydiffnorm;
    VecCopy(Y, Ylast);
    VecScale(Y, alpha);
    if (precheck->monitor) {
      PetscViewerASCIIPrintf(
          precheck->monitor,
          "Angle %g degrees less than threshold %g, corrected step by alpha=%g\n",
          (double)(theta * 180. / PETSC_PI),
          (double)precheck->angle,
          (double)alpha);
    }
  } else {
    VecCopy(Y, Ylast);
    *changed = PETSC_FALSE;
    if (precheck->monitor) {
      PetscViewerASCIIPrintf(
          precheck->monitor,
          "Angle %g degrees exceeds threshold %g, no correction applied\n",
          (double)(theta * 180. / PETSC_PI),
          (double)precheck->angle);
    }
  }
  return 0;
}
PetscErrorCode PreCheckDestroy(PreCheck* precheck) {
  if (!*precheck)
    return 0;
  VecDestroy(&(*precheck)->Ylast);
  PetscViewerDestroy(&(*precheck)->monitor);
  PetscFree(*precheck);
  return 0;
}
PetscErrorCode PreCheckCreate(MPI_Comm comm, PreCheck* precheck) {
  PetscNew(precheck);
  (*precheck)->comm  = comm;
  (*precheck)->angle = 10.; /* only active if angle is less than 10 degrees */
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
                                uy_S = dhy * (u - u_S), e_E = eta(user, i, j, ux_E, uy_E),
                                e_W = eta(user, i, j, ux_W, uy_W), e_N = eta(user, i, j, ux_N, uy_N),
                                e_S = eta(user, i, j, ux_S, uy_S), de_E = deta(user, xx, yy, ux_E, uy_E),
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
