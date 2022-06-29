#include "transient_groundwater.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscts.h>

/*
  User-defined routines and data structures
*/
typedef struct {
  PetscScalar u, v, omega, temp;
} Field;

namespace FanDarcyGroundwater {
PetscErrorCode FormIFunctionLocal(DMDALocalInfo*, PetscReal, Field**, Field**, Field**, void*);
PetscErrorCode FormInitialSolution(TS, Vec, AppCtx*);

int update(Parameters& params, ArrayPack& arp, AppCtx& user_context, DMDA_Array_Pack& dmdapack) {
  PetscInt mx, my, steps;
  TS ts;
  Vec X;
  PetscReal ftime;
  TSConvergedReason reason;

  TSCreate(PETSC_COMM_WORLD, &ts);
  DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      DMDA_STENCIL_STAR,
      4,
      4,
      PETSC_DECIDE,
      PETSC_DECIDE,
      4,
      1,
      0,
      0,
      &user_context.da);
  DMSetFromOptions(user_context.da);
  DMSetUp(user_context.da);
  TSSetDM(ts, user_context.da);

  DMDAGetInfo(
      user_context.da,
      0,
      &mx,
      &my,
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
     Problem parameters (velocity of lid, prandtl, and grashof numbers)
  */
  user_context.lidvelocity = 1.0 / (mx * my);
  user_context.prandtl     = 1.0;
  user_context.grashof     = 1.0;
  user_context.parabolic   = PETSC_FALSE;
  user_context.cfl_initial = 50.;

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Driven cavity/natural convection options", "");
  PetscOptionsReal(
      "-lidvelocity",
      "Lid velocity, related to Reynolds number",
      "",
      user_context.lidvelocity,
      &user_context.lidvelocity,
      NULL);
  PetscOptionsReal(
      "-prandtl", "Ratio of viscous to thermal diffusivity", "", user_context.prandtl, &user_context.prandtl, NULL);
  PetscOptionsReal(
      "-grashof", "Ratio of bouyant to viscous forces", "", user_context.grashof, &user_context.grashof, NULL);
  PetscOptionsBool(
      "-parabolic",
      "Relax incompressibility to make the system parabolic instead of differential-algebraic",
      "",
      user_context.parabolic,
      &user_context.parabolic,
      NULL);
  PetscOptionsReal(
      "-cfl_initial",
      "Advective CFL for the first time step",
      "",
      user_context.cfl_initial,
      &user_context.cfl_initial,
      NULL);
  PetscOptionsEnd();

  DMDASetFieldName(user_context.da, 0, "x-velocity");
  DMDASetFieldName(user_context.da, 1, "y-velocity");
  DMDASetFieldName(user_context.da, 2, "Omega");
  DMDASetFieldName(user_context.da, 3, "temperature");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create user context, set problem data, create vector data structures.
     Also, compute the initial guess.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create time integration context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMSetApplicationContext(user_context.da, &user_context);
  DMDATSSetIFunctionLocal(user_context.da, INSERT_VALUES, (DMDATSIFunctionLocal)FormIFunctionLocal, &user_context);
  TSSetMaxSteps(ts, 10000);
  TSSetMaxTime(ts, 1.0);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);
  TSSetTimeStep(ts, 0.1);
  TSSetFromOptions(ts);

  PetscPrintf(
      PETSC_COMM_WORLD,
      "%" PetscInt_FMT "x%" PetscInt_FMT " grid, lid velocity = %g, prandtl # = %g, grashof # = %g\n",
      mx,
      my,
      (double)user_context.lidvelocity,
      (double)user_context.prandtl,
      (double)user_context.grashof);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  DMCreateGlobalVector(user_context.da, &X);
  FormInitialSolution(ts, X, &user_context);

  TSSolve(ts, X);
  TSGetSolveTime(ts, &ftime);
  TSGetStepNumber(ts, &steps);
  TSGetConvergedReason(ts, &reason);

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
  VecDestroy(&X);
  TSDestroy(&ts);

  return 0;
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
PetscErrorCode FormInitialSolution(TS ts, Vec X, AppCtx* user_context) {
  DM da;
  PetscInt i, j, mx, xs, ys, xm, ym;
  PetscReal grashof, dx;
  PetscScalar **x, **my_topo;

  grashof = user_context->grashof;
  TSGetDM(ts, &da);
  DMDAGetInfo(da, 0, &mx, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  dx = 1.0 / (mx - 1);

  /*
     Get local grid boundaries (for 2-dimensional DMDA):
       xs, ys   - starting grid indices (no ghost points)
       xm, ym   - widths of local grid (no ghost points)
  */
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);

  /*
     Get a pointer to vector data.
       - For default PETSc vectors, VecGetArray() returns a pointer to
         the data array.  Otherwise, the routine is implementation dependent.
       - You MUST call VecRestoreArray() when you no longer need access to
         the array.
  */
  DMDAVecGetArray(da, X, &x);
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));  // TODO: this is causing current issue

  /*
     Compute initial guess over the locally owned part of the grid
     Initial condition is motionless fluid and equilibrium temperature
  */
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      x[j][i] = 0.0;
    }
  }

  /*
     Restore vector
  */
  DMDAVecRestoreArray(da, X, &x);
  PetscCall(DMDAVecRestoreArray(da, user_context->topo_vec, &my_topo));

  return 0;
}

PetscErrorCode FormIFunctionLocal(DMDALocalInfo* info, PetscReal ptime, Field** x, Field** xdot, Field** f, void* ptr) {
  AppCtx* user_context = (AppCtx*)ptr;
  PetscInt xints, xinte, yints, yinte, i, j;
  PetscReal hx, hy, dhx, dhy, hxdhy, hydhx;
  PetscReal grashof, prandtl, lid;
  PetscScalar u, udot, uxx, uyy, vx, vy, avx, avy, vxp, vxm, vyp, vym;

  grashof = user_context->grashof;
  prandtl = user_context->prandtl;
  lid     = user_context->lidvelocity;

  /*
     Define mesh intervals ratios for uniform grid.

   Note: FD formulae below are normalized by multiplying through by
   local volume element (i.e. hx*hy) to obtain coefficients O(1) in two dimensions.

*/
  dhx   = (PetscReal)(info->mx - 1);
  dhy   = (PetscReal)(info->my - 1);
  hx    = 1.0 / dhx;
  hy    = 1.0 / dhy;
  hxdhy = hx * dhy;
  hydhx = hy * dhx;

  xints = info->xs;
  xinte = info->xs + info->xm;
  yints = info->ys;
  yinte = info->ys + info->ym;

  /* Test whether we are on the bottom edge of the global array */
  if (yints == 0) {
    j     = 0;
    yints = yints + 1;
    /* bottom edge */
    for (i = info->xs; i < info->xs + info->xm; i++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega + (x[j + 1][i].u - x[j][i].u) * dhy;
      f[j][i].temp  = x[j][i].temp - x[j + 1][i].temp;
    }
  }

  /* Test whether we are on the top edge of the global array */
  if (yinte == info->my) {
    j     = info->my - 1;
    yinte = yinte - 1;
    /* top edge */
    for (i = info->xs; i < info->xs + info->xm; i++) {
      f[j][i].u     = x[j][i].u - lid;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega + (x[j][i].u - x[j - 1][i].u) * dhy;
      f[j][i].temp  = x[j][i].temp - x[j - 1][i].temp;
    }
  }

  /* Test whether we are on the left edge of the global array */
  if (xints == 0) {
    i     = 0;
    xints = xints + 1;
    /* left edge */
    for (j = info->ys; j < info->ys + info->ym; j++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega - (x[j][i + 1].v - x[j][i].v) * dhx;
      f[j][i].temp  = x[j][i].temp;
    }
  }

  /* Test whether we are on the right edge of the global array */
  if (xinte == info->mx) {
    i     = info->mx - 1;
    xinte = xinte - 1;
    /* right edge */
    for (j = info->ys; j < info->ys + info->ym; j++) {
      f[j][i].u     = x[j][i].u;
      f[j][i].v     = x[j][i].v;
      f[j][i].omega = x[j][i].omega - (x[j][i].v - x[j][i - 1].v) * dhx;
      f[j][i].temp  = x[j][i].temp - (PetscReal)(grashof > 0);
    }
  }

  /* Compute over the interior points */
  for (j = yints; j < yinte; j++) {
    for (i = xints; i < xinte; i++) {
      /*
        convective coefficients for upwinding
      */
      vx  = x[j][i].u;
      avx = PetscAbsScalar(vx);
      vxp = .5 * (vx + avx);
      vxm = .5 * (vx - avx);
      vy  = x[j][i].v;
      avy = PetscAbsScalar(vy);
      vyp = .5 * (vy + avy);
      vym = .5 * (vy - avy);

      /* U velocity */
      u         = x[j][i].u;
      udot      = user_context->parabolic ? xdot[j][i].u : 0.;
      uxx       = (2.0 * u - x[j][i - 1].u - x[j][i + 1].u) * hydhx;
      uyy       = (2.0 * u - x[j - 1][i].u - x[j + 1][i].u) * hxdhy;
      f[j][i].u = udot + uxx + uyy - .5 * (x[j + 1][i].omega - x[j - 1][i].omega) * hx;

      /* V velocity */
      u         = x[j][i].v;
      udot      = user_context->parabolic ? xdot[j][i].v : 0.;
      uxx       = (2.0 * u - x[j][i - 1].v - x[j][i + 1].v) * hydhx;
      uyy       = (2.0 * u - x[j - 1][i].v - x[j + 1][i].v) * hxdhy;
      f[j][i].v = udot + uxx + uyy + .5 * (x[j][i + 1].omega - x[j][i - 1].omega) * hy;

      /* Omega */
      u   = x[j][i].omega;
      uxx = (2.0 * u - x[j][i - 1].omega - x[j][i + 1].omega) * hydhx;
      uyy = (2.0 * u - x[j - 1][i].omega - x[j + 1][i].omega) * hxdhy;
      f[j][i].omega =
          (xdot[j][i].omega + uxx + uyy + (vxp * (u - x[j][i - 1].omega) + vxm * (x[j][i + 1].omega - u)) * hy +
           (vyp * (u - x[j - 1][i].omega) + vym * (x[j + 1][i].omega - u)) * hx -
           .5 * grashof * (x[j][i + 1].temp - x[j][i - 1].temp) * hy);

      /* Temperature */
      u   = x[j][i].temp;
      uxx = (2.0 * u - x[j][i - 1].temp - x[j][i + 1].temp) * hydhx;
      uyy = (2.0 * u - x[j - 1][i].temp - x[j + 1][i].temp) * hxdhy;
      f[j][i].temp =
          (xdot[j][i].temp + uxx + uyy +
           prandtl * ((vxp * (u - x[j][i - 1].temp) + vxm * (x[j][i + 1].temp - u)) * hy +
                      (vyp * (u - x[j - 1][i].temp) + vym * (x[j + 1][i].temp - u)) * hx));
    }
  }
  std::cout << "x in the func " << x[1][1].u << " f " << f[1][1].u << std::endl;

  /*
     Flop count (multiply-adds are counted as 2 operations)
  */
  PetscLogFlops(84.0 * info->ym * info->xm);
  return 0;
}

}  // namespace FanDarcyGroundwater
