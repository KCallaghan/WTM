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

PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);
PetscErrorCode FormFunction(SNES, Vec, Vec, void*);
PetscErrorCode FormInitialGuess(Vec, Parameters& params, ArrayPack& arp);
PetscErrorCode Monitor(SNES, PetscInt, PetscReal, void*);
PetscErrorCode MatrixFreePreconditioner(PC, Vec, Vec);

/*
   User-defined application context
*/
typedef struct {
  DM da; /* distributed array */
  Vec F; /* right-hand-side of PDE */
  Vec ctx_QE, ctx_QN, ctx_QW, ctx_QS, ctx_head, ctx_storativity;
  PetscMPIInt rank; /* rank of processor */
  PetscMPIInt size; /* size of communicator */
  PetscReal h;      /* mesh spacing */
  PetscBool sjerr;  /* if or not to test jacobian domain error */
  PetscReal timestep;
  PetscReal ncells_x;
} ApplicationCtx;

/*
   User-defined context for monitoring
*/
typedef struct {
  PetscViewer viewer;
} MonitorCtx;

/*
   User-defined context for checking candidate iterates that are
   determined by line search methods
*/
typedef struct {
  Vec last_step;       /* previous iterate */
  PetscReal tolerance; /* tolerance for changes between successive iterates */
  ApplicationCtx* user;
} StepCheckCtx;

typedef struct {
  PetscInt its0; /* num of prevous outer KSP iterations */
} SetSubKSPCtx;

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

#pragma omp parallel for default(none) shared(arp, params) collapse(2)
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
    }
  }
}

int update(Parameters& params, ArrayPack& arp) {
  SNES snes; /* SNES context */
  KSP ksp;
  PC pc;
  SNESLineSearch linesearch;                                                 /* SNESLineSearch context */
  Mat J;                                                                     /* Jacobian matrix */
  ApplicationCtx ctx;                                                        /* user-defined context */
  Vec x, r, U, F, ctx_QE, ctx_QN, ctx_QW, ctx_QS, ctx_head, ctx_storativity; /* vectors */
  MonitorCtx monP;                                                           /* monitoring context */
  StepCheckCtx checkP;                                                       /* step-checking context */
  SetSubKSPCtx checkP1;
  PetscBool pre_check, post_check, post_setsubksp; /* flag indicating whether we're checking candidate iterates */
  PetscScalar xp, *FF, *UU, none = -1.0, *QE, *QW, *QN, *QS, *storativity, *head, *my_x;
  PetscInt its, N = params.ncells_x * params.ncells_y, i, maxit, maxf, xs,
                xm;  // N seems like it might be the total number of cells in the vectors
  PetscReal abstol, rtol, stol, norm;
  PetscBool flg, viewinitial = PETSC_FALSE;

  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &ctx.rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &ctx.size));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &N, NULL));
  ctx.h        = 1.0 / (N - 1);
  ctx.sjerr    = PETSC_FALSE;
  ctx.timestep = params.deltat;
  ctx.ncells_x = params.ncells_x;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-test_jacobian_domain_error", &ctx.sjerr, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-view_initial", &viewinitial, NULL));

  // set starting values and arrays for the groundwater calculation
  set_starting_values(params, arp);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create vector data structures; set function evaluation routine
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create distributed array (DMDA) to manage parallel grid and vectors
  */
  PetscCall(DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, N, 1, 1, NULL, &ctx.da));
  PetscCall(DMSetFromOptions(ctx.da));
  PetscCall(DMSetUp(ctx.da));

  /*
     Extract global and local vectors from DMDA; then duplicate for remaining
     vectors that are the same types
  */
  PetscCall(DMCreateGlobalVector(ctx.da, &x));
  PetscCall(VecDuplicate(x, &r));
  PetscCall(VecDuplicate(x, &F));
  ctx.F = F;
  PetscCall(VecDuplicate(x, &U));
  PetscCall(VecDuplicate(x, &ctx_QE));
  ctx.ctx_QE = ctx_QE;
  PetscCall(VecDuplicate(x, &ctx_QW));
  ctx.ctx_QW = ctx_QW;
  PetscCall(VecDuplicate(x, &ctx_QN));
  ctx.ctx_QN = ctx_QN;
  PetscCall(VecDuplicate(x, &ctx_QS));
  ctx.ctx_QS = ctx_QS;
  PetscCall(VecDuplicate(x, &ctx_storativity));
  ctx.ctx_storativity = ctx_storativity;
  PetscCall(VecDuplicate(x, &ctx_head));
  ctx.ctx_head = ctx_head;

  PetscCall(DMDAVecGetArray(ctx.da, ctx_QE, &QE));
  PetscCall(DMDAVecGetArray(ctx.da, ctx_QW, &QW));
  PetscCall(DMDAVecGetArray(ctx.da, ctx_QN, &QN));
  PetscCall(DMDAVecGetArray(ctx.da, ctx_QS, &QS));
  PetscCall(DMDAVecGetArray(ctx.da, ctx_storativity, &storativity));
  PetscCall(DMDAVecGetArray(ctx.da, ctx_head, &head));

  // TODO: check how to set up the edge cells correctly
  for (int count_y = 1; count_y < params.ncells_y - 1; count_y++)
    for (int count_x = 1; count_x < params.ncells_x - 1; count_x++) {
      double edge_T = 2 / (1 / arp.transmissivity(count_x, count_y) + 1 / arp.transmissivity(count_x, count_y + 1));
      int count_i   = arp.wtd.xyToI(count_x, count_y);
      QE[count_i]   = edge_T *
                    ((arp.wtd(count_x, count_y) + arp.topo(count_x, count_y)) -
                     (arp.wtd(count_x, count_y + 1) + arp.topo(count_x, count_y + 1))) /
                    (arp.cellsize_e_w_metres[count_y] * arp.cellsize_e_w_metres[count_y]);
      edge_T      = 2 / (1 / arp.transmissivity(count_x, count_y) + 1 / arp.transmissivity(count_x, count_y - 1));
      QW[count_i] = edge_T *
                    ((arp.wtd(count_x, count_y) + arp.topo(count_x, count_y)) -
                     (arp.wtd(count_x, count_y - 1) + arp.topo(count_x, count_y - 1))) /
                    (arp.cellsize_e_w_metres[count_y] * arp.cellsize_e_w_metres[count_y]);
      edge_T      = 2 / (1 / arp.transmissivity(count_x, count_y) + 1 / arp.transmissivity(count_x - 1, count_y));
      QN[count_i] = edge_T *
                    ((arp.wtd(count_x, count_y) + arp.topo(count_x, count_y)) -
                     (arp.wtd(count_x - 1, count_y) + arp.topo(count_x - 1, count_y))) /
                    (params.cellsize_n_s_metres * params.cellsize_n_s_metres);
      edge_T      = 2 / (1 / arp.transmissivity(count_x, count_y) + 1 / arp.transmissivity(count_x + 1, count_y));
      QS[count_i] = edge_T *
                    ((arp.wtd(count_x, count_y) + arp.topo(count_x, count_y)) -
                     (arp.wtd(count_x + 1, count_y) + arp.topo(count_x + 1, count_y))) /
                    (params.cellsize_n_s_metres * params.cellsize_n_s_metres);

      storativity[count_i] = arp.effective_storativity(count_x, count_y);
      head[count_i]        = arp.wtd(count_x, count_y) + arp.topo(count_x, count_y);
    }

  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_QE, &QE));
  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_QW, &QW));
  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_QN, &QN));
  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_QS, &QS));
  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_storativity, &storativity));
  PetscCall(DMDAVecRestoreArray(ctx.da, ctx_head, &head));

  /*
     Set function evaluation routine and vector.  Whenever the nonlinear
     solver needs to compute the nonlinear function, it will call this
     routine.
      - Note that the final routine argument is the user-defined
        context that provides application-specific data for the
        function evaluation routine.
  */

  PetscCall(SNESSetFunction(snes, r, FormFunction, &ctx));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; set Jacobian evaluation routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(MatCreate(PETSC_COMM_WORLD, &J));
  PetscCall(MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatSetFromOptions(J));
  PetscCall(MatSeqAIJSetPreallocation(J, 5, NULL));
  PetscCall(MatMPIAIJSetPreallocation(J, 3, NULL, 2, NULL));

  /*
     Set Jacobian matrix data structure and default Jacobian evaluation
     routine.  Whenever the nonlinear solver needs to compute the
     Jacobian matrix, it will call this routine.
      - Note that the final routine argument is the user-defined
        context that provides application-specific data for the
        Jacobian evaluation routine.
  */
  PetscCall(SNESSetJacobian(snes, J, J, FormJacobian, &ctx));

  /*
     Optionally allow user-provided preconditioner
   */

  PetscCall(PetscOptionsHasName(NULL, NULL, "-user_precond", &flg));
  if (flg) {
    KSP ksp;
    PC pc;
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCSHELL));
    PetscCall(PCShellSetApply(pc, MatrixFreePreconditioner));
  }

  //     flg  = PETSC_FALSE;
  //   PetscOptionsGetBool(NULL,NULL,"-snes_mf",&flg,NULL);
  //
  // // PetscCall(PetscOptionsHasName(NULL, NULL, "-user_precond", &flg));
  //  if (flg) {
  //     SNESGetKSP(snes,&ksp);
  //     KSPGetPC(ksp,&pc);
  //     PetscOptionsHasName(NULL,NULL,"-user_precond",&flg);
  //     if (flg) { /* user-defined precond */
  //       PCSetType(pc,PCSHELL);
  //       PCShellSetApply(pc,MatrixFreePreconditioner);
  //     } else {PCSetType(pc,PCNONE);}
  //   }
  //   SNESSetFromOptions(snes);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Set an optional user-defined monitoring routine
  */
  PetscCall(PetscViewerDrawOpen(PETSC_COMM_WORLD, 0, 0, 0, 0, 400, 400, &monP.viewer));
  PetscCall(SNESMonitorSet(snes, Monitor, &monP, 0));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize application:
     Store right-hand-side of PDE and exact solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Get local grid boundaries (for 1-dimensional DMDA):
       xs, xm - starting grid index, width of local grid (no ghost points)
  */
  PetscCall(DMDAGetCorners(ctx.da, &xs, NULL, NULL, &xm, NULL, NULL));

  /*
     Get pointers to vector data
  */
  PetscCall(DMDAVecGetArray(ctx.da, F, &FF));
  PetscCall(DMDAVecGetArray(ctx.da, U, &UU));

  /*
     Compute local vector entries
  */
  xp = ctx.h * xs;
  for (i = xs; i < xs + xm; i++) {
    FF[i] = 6.0 * xp + PetscPowScalar(xp + 1.e-12, 6.0); /* +1.e-12 is to prevent 0^6 */
    UU[i] = xp * xp * xp;
    xp += ctx.h;
  }

  /*
     Restore vectors
  */
  PetscCall(DMDAVecRestoreArray(ctx.da, F, &FF));
  PetscCall(DMDAVecRestoreArray(ctx.da, U, &UU));
  if (viewinitial) {
    PetscCall(VecView(U, 0));
    PetscCall(VecView(F, 0));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Evaluate initial guess; then solve nonlinear system
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Note: The user should initialize the vector, x, with the initial guess
     for the nonlinear solver prior to calling SNESSolve().  In particular,
     to employ an initial guess of zero, the user should explicitly set
     this vector to zero by calling VecSet().
  */
  PetscCall(FormInitialGuess(x, params, arp));
  PetscCall(SNESSolve(snes, NULL, x));
  PetscCall(SNESGetIterationNumber(snes, &its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Number of SNES iterations = %" PetscInt_FMT "\n", its));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check solution and clean up
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  PetscCall(DMDAVecGetArray(ctx.da, x, &my_x));

  for (int count_y = 0; count_y < params.ncells_y; count_y++)
    for (int count_x = 0; count_x < params.ncells_x; count_x++) {
      int count_i               = arp.wtd.xyToI(count_x, count_y);
      double my_head            = my_x[count_i];
      arp.wtd(count_x, count_y) = my_head - arp.topo(count_x, count_y);
    }
  PetscCall(DMDAVecRestoreArray(ctx.da, x, &my_x));

  PetscCall(PetscViewerDestroy(&monP.viewer));
  if (post_check)
    PetscCall(VecDestroy(&checkP.last_step));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&U));
  PetscCall(VecDestroy(&F));
  PetscCall(MatDestroy(&J));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&ctx.da));

  for (int count_y = 0; count_y < params.ncells_y; count_y++) {
    for (int count_x = 0; count_x < params.ncells_x; count_x++) {
      if (arp.land_mask(count_x, count_y) != 0) {
        continue;
      }
      // count up any water lost to the ocean so that we can compute a water balance
      if (arp.wtd(count_x, count_y) > 0) {
        arp.total_loss_to_ocean += arp.wtd(count_x, count_y) * arp.cell_area[count_y];
      } else {
        arp.total_loss_to_ocean += arp.wtd(count_x, count_y) * arp.cell_area[count_y] * arp.porosity(count_x, count_y);
      }
      //   arp.wtd(count_x, count_y) = 0.;  // reset ocean water tables to 0.
    }
  }

  return 0;
}

/* ------------------------------------------------------------------- */
/*
   FormInitialGuess - Computes initial guess.

   Input/Output Parameter:
.  x - the solution vector
*/
PetscErrorCode FormInitialGuess(Vec x, Parameters& params, ArrayPack& arp) {
  // I think this is where we put in the wtd + topo as the starting guess for the solver
  PetscFunctionBeginUser;

  for (int count_y = 0; count_y < params.ncells_y; count_y++)
    for (int count_x = 0; count_x < params.ncells_x; count_x++) {
      int count_i    = arp.wtd.xyToI(count_x, count_y);
      double my_head = arp.wtd(count_x, count_y) + arp.topo(count_x, count_y);
      VecSetValue(x, count_i, my_head, INSERT_VALUES);
    }

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   FormFunction - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  ctx - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  f - function vector

   Note:
   The user-defined context can contain any application-specific
   data needed for the function evaluation.
*/
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx) {
  ApplicationCtx* user = (ApplicationCtx*)ctx;
  DM da                = user->da;
  double timestep      = user->timestep;
  PetscScalar *xx, *ff, *FF, *QE, *QW, *QN, *QS, *storativity, *head;  //,d;
  PetscInt i, M, xs, xm;
  Vec xlocal;

  PetscFunctionBeginUser;
  PetscCall(DMGetLocalVector(da, &xlocal));
  /*
     Scatter ghost points to local vector, using the 2-step process
        DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
     By placing code between these two statements, computations can
     be done while messages are in transition.
  */
  PetscCall(DMGlobalToLocalBegin(da, x, INSERT_VALUES, xlocal));
  PetscCall(DMGlobalToLocalEnd(da, x, INSERT_VALUES, xlocal));

  /*
     Get pointers to vector data.
       - The vector xlocal includes ghost point; the vectors x and f do
         NOT include ghost points.
       - Using DMDAVecGetArray() allows accessing the values using global ordering
  */
  PetscCall(DMDAVecGetArray(da, xlocal, &xx));
  PetscCall(DMDAVecGetArray(da, f, &ff));
  PetscCall(DMDAVecGetArray(da, user->F, &FF));
  PetscCall(DMDAVecGetArray(da, user->ctx_QE, &QE));
  PetscCall(DMDAVecGetArray(da, user->ctx_QW, &QW));
  PetscCall(DMDAVecGetArray(da, user->ctx_QN, &QN));
  PetscCall(DMDAVecGetArray(da, user->ctx_QS, &QS));
  PetscCall(DMDAVecGetArray(da, user->ctx_storativity, &storativity));
  PetscCall(DMDAVecGetArray(da, user->ctx_head, &head));

  /*
     Get local grid boundaries (for 1-dimensional DMDA):
       xs, xm  - starting grid index, width of local grid (no ghost points)
  */
  PetscCall(DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL));
  PetscCall(DMDAGetInfo(da, NULL, &M, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));

  /*
     Set function values for boundary points; define local interior grid point range:
        xsi - starting interior grid index
        xei - ending interior grid index
  */
  if (xs == 0) { /* left boundary */
    ff[0] = xx[0];
    xs++;
    xm--;
  }
  if (xs + xm == M) { /* right boundary */
    ff[xs + xm - 1] = xx[xs + xm - 1] - 1.0;
    xm--;
  }

  /*
     Compute function over locally owned part of the grid (interior points only)
  */

  // d = 1.0/(user->h*user->h);
  for (i = xs; i < xs + xm; i++)
    ff[i] = storativity[i] * (xx[i] - head[i]) / timestep - QS[i] + QN[i] - QE[i] + QW[i];

  /*
     Restore vectors
  */
  PetscCall(DMDAVecRestoreArray(da, xlocal, &xx));
  PetscCall(DMDAVecRestoreArray(da, f, &ff));
  PetscCall(DMDAVecRestoreArray(da, user->F, &FF));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_QE, &QE));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_QW, &QW));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_QN, &QN));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_QS, &QS));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_storativity, &storativity));
  PetscCall(DMDAVecRestoreArray(da, user->ctx_head, &head));

  PetscCall(DMRestoreLocalVector(da, &xlocal));
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   FormJacobian - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void* ctx) {
  ApplicationCtx* user = (ApplicationCtx*)ctx;
  PetscScalar *xx, d, A[3];
  PetscInt i, j[3], M, xs, xm;
  DM da = user->da;

  PetscFunctionBeginUser;
  if (user->sjerr) {
    PetscCall(SNESSetJacobianDomainError(snes));
    PetscFunctionReturn(0);
  }
  /*
     Get pointer to vector data
  */
  PetscCall(DMDAVecGetArrayRead(da, x, &xx));
  PetscCall(DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL));

  /*
    Get range of locally owned matrix
  */
  PetscCall(DMDAGetInfo(da, NULL, &M, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));

  /*
     Determine starting and ending local indices for interior grid points.
     Set Jacobian entries for boundary points.
  */

  if (xs == 0) { /* left boundary */
    i    = 0;
    A[0] = 1.0;

    PetscCall(MatSetValues(jac, 1, &i, 1, &i, A, INSERT_VALUES));
    xs++;
    xm--;
  }
  if (xs + xm == M) { /* right boundary */
    i    = M - 1;
    A[0] = 1.0;
    PetscCall(MatSetValues(jac, 1, &i, 1, &i, A, INSERT_VALUES));
    xm--;
  }

  /*
     Interior grid points
      - Note that in this case we set all elements for a particular
        row at once.
  */
  d = 1.0 / (user->h * user->h);
  for (i = xs; i < xs + xm; i++) {
    j[0] = i - 1;
    j[1] = i;
    j[2] = i + 1;
    A[0] = A[2] = d;
    A[1]        = -2.0 * d + 2.0 * xx[i];
    PetscCall(MatSetValues(jac, 1, &i, 3, j, A, INSERT_VALUES));
  }

  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.

     Also, restore vector.
  */

  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(DMDAVecRestoreArrayRead(da, x, &xx));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   Monitor - Optional user-defined monitoring routine that views the
   current iterate with an x-window plot. Set by SNESMonitorSet().

   Input Parameters:
   snes - the SNES context
   its - iteration number
   norm - 2-norm function value (may be estimated)
   ctx - optional user-defined context for private data for the
         monitor routine, as set by SNESMonitorSet()

   Note:
   See the manpage for PetscViewerDrawOpen() for useful runtime options,
   such as -nox to deactivate all x-window output.
 */
PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal fnorm, void* ctx) {
  MonitorCtx* monP = (MonitorCtx*)ctx;
  Vec x;

  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iter = %" PetscInt_FMT ",SNES Function norm %g\n", its, (double)fnorm));
  PetscCall(SNESGetSolution(snes, &x));
  PetscCall(VecView(x, monP->viewer));
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/*
   MatrixFreePreconditioner - This routine demonstrates the use of a
   user-provided preconditioner.  This code implements just the null
   preconditioner, which of course is not recommended for general use.

   Input Parameters:
+  pc - preconditioner
-  x - input vector

   Output Parameter:
.  y - preconditioned vector
*/
PetscErrorCode MatrixFreePreconditioner(PC pc, Vec x, Vec y) {
  PetscCall(VecCopy(x, y));
  return 0;
}

}  // namespace FanDarcyGroundwater
