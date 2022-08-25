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

#include <Eigen/Core>
#include <Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;  // declares a row-major sparse matrix type of double

constexpr double seconds_in_a_year = 31536000.;

std::string my_time() {
  const auto now       = std::chrono::system_clock::now();
  const auto in_time_t = std::chrono::system_clock::to_time_t(now);

  std::stringstream ss;
  ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %H:%M:%S");
  return ss.str();
}

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
      } else {
        double rech_count = arp.rech(x, y);
        if (arp.wtd(x, y) >= 0 && arp.wtd(x, y) + arp.rech(x, y) < 0)
          rech_count = -arp.wtd(x, y);

        arp.total_added_recharge += rech_count * arp.cell_area[y];
      }
    }
  }

#pragma omp parallel for default(none) shared(arp, params) collapse(2)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask(x, y) == 0.f || arp.wtd(x, y) >= 0) {
        arp.effective_storativity(x, y) = 1.;
      } else {
        arp.effective_storativity(x, y) = arp.porosity(x, y);
      }

      if (arp.land_mask(x, y) == 0.f) {
        // in the ocean, we set several arrays to default values
        arp.transmissivity(x, y)  = ocean_T;
        arp.wtd_T(x, y)           = 0.;
        arp.original_wtd(x, y)    = 0.;
        arp.wtd_T_iteration(x, y) = 0.;
      } else {
        arp.transmissivity(x, y) = depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y));
        arp.original_wtd(x, y)   = arp.wtd(x, y);
        // Apply the first half of the recharge to the water-table depth grid (wtd)
        // Its clone (wtd_T) is used and updated in the Picard iteration
        // use regular porosity for adding recharge since this checks
        // for underground space within add_recharge.
        arp.wtd(x, y) += add_recharge(arp.rech(x, y) / 5., arp.wtd(x, y), arp.porosity(x, y));

        arp.wtd_T(x, y)           = arp.wtd(x, y);
        arp.wtd_T_iteration(x, y) = arp.wtd_T(x, y);
      }
      // set the scalar arrays for x and y directions
      arp.scalar_array_y(x, y) = params.deltat / (5. * arp.effective_storativity(x, y) * arp.cellsize_e_w_metres[y] *
                                                  arp.cellsize_e_w_metres[y]);
      params.x_partial         = params.deltat / (5. * params.cellsize_n_s_metres * params.cellsize_n_s_metres);
      arp.scalar_array_x(x, y) = params.x_partial / arp.effective_storativity(x, y);
    }
  }
}

// The midpoint method consists of two main steps.
// In the first step, we compute water tables at the time that is *half* of the full time-step.
// To do so, we use an implicit backward-difference Euler method.
// Using these half-time water tables, we compute the new transmissivity values
// That will be used in the second step of the midpoint method below.
void first_half(const Parameters& params, ArrayPack& arp) {
  Eigen::VectorXd b(params.ncells_x * params.ncells_y);
  Eigen::VectorXd vec_x(params.ncells_x * params.ncells_y);
  Eigen::VectorXd guess(params.ncells_x * params.ncells_y);
  SpMat A(params.ncells_x * params.ncells_y, params.ncells_x * params.ncells_y);
  // reserve the needed space in the A matrix. We know there is a maximum of 5 items per matrix row or column.
  A.reserve(Eigen::VectorXi::Constant(params.ncells_x * params.ncells_y, 5));

// We need to solve the vector-matrix equation Ax=b.
// b consists of the current head values (i.e. water table depth + topography)
// We also populate a 'guess', which consists of a water table (wtd_T) that in later iterations
// has already been modified for changing transmissivity closer to the final answer.
#pragma omp parallel for default(none) shared(arp, b, guess, params) collapse(2)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      arp.wtd_T_iteration(x, y) = arp.wtd_T(x, y);
      b(arp.wtd_T.xyToI(x, y)) =
          arp.wtd(x, y) + arp.topo(x, y);  // wtd is 0 in ocean cells and topo is 0 in ocean cells, so no need to
                                           // differentiate between ocean vs land.
      guess(arp.wtd_T.xyToI(x, y)) = arp.wtd_T(x, y) + arp.topo(x, y);
    }
  }

  //  HALFWAY SOLVE
  //  populate the A matrix.

  const auto construct_e_w_diagonal_one = [&](const int x, const int y, const int dy) {
    return -arp.scalar_array_y(x, y) * ((arp.transmissivity(x, y + dy) + arp.transmissivity(x, y)) / 2);
  };

  const auto construct_n_s_diagonal_one = [&](const int x, const int y, const int dx) {
    return -arp.scalar_array_x(x, y) * ((arp.transmissivity(x + dx, y) + arp.transmissivity(x, y)) / 2);
  };

  const auto construct_major_diagonal_one =
      [&](const int x, const int y, const int dx1, const int dx2, const int dy1, const int dy2) {
        const auto x_term = arp.scalar_array_x(x, y) * (arp.transmissivity(x + dx1, y) / 2 + arp.transmissivity(x, y) +
                                                        arp.transmissivity(x + dx2, y) / 2);
        const auto y_term = arp.scalar_array_y(x, y) * (arp.transmissivity(x, y + dy1) / 2 + arp.transmissivity(x, y) +
                                                        arp.transmissivity(x, y + dy2) / 2);
        return x_term + y_term + 1;
      };

#pragma omp parallel for collapse(2) default(none) \
    shared(arp, params, construct_major_diagonal_one, construct_e_w_diagonal_one, construct_n_s_diagonal_one, A)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      // The row and column that the current cell will be stored in in matrix A.
      // This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
      // All of the N,E,S,W directions should be in the same row, but the column will differ.
      const auto main_loc = arp.wtd_T.xyToI(x, y);

      if (x != 0) {
        // Do the North diagonal. Offset by -1. When x == 0, there is no north diagonal.
        A.insert(main_loc, main_loc - 1) = construct_n_s_diagonal_one(x, y, -1);
      }

      if (y != 0) {
        // Next is the West diagonal. Opposite of the East. Located at (i,j-params.ncells_x). When y == 0, there is no
        // west diagonal.
        A.insert(main_loc, main_loc - params.ncells_x) = construct_e_w_diagonal_one(x, y, -1);
      }

      // major diagonal:
      A.insert(main_loc, main_loc) = construct_major_diagonal_one(
          x,
          y,
          (x == 0) ? 0 : -1,
          (x == params.ncells_x - 1) ? 0 : 1,
          (y == 0) ? 0 : -1,
          (y == params.ncells_y - 1) ? 0 : 1);

      if (y != params.ncells_y - 1) {
        // Now do the East diagonal. The East location is at (i,j+params.ncells_x). When y == params.ncells_y -1,
        // there is no east diagonal.
        A.insert(main_loc, main_loc + params.ncells_x) = construct_e_w_diagonal_one(x, y, 1);
      }

      if (x != params.ncells_x - 1) {
        // Do the South diagonal, offset by +1. When x == params.ncells_x, there is no south diagonal.
        A.insert(main_loc, main_loc + 1) = construct_n_s_diagonal_one(x, y, 1);
      }
    }
  }

  A.makeCompressed();

  // Biconjugate gradient solver with guess
  Eigen::BiCGSTAB<SpMat> solver;
  solver.setTolerance(params.solver_tolerance_value);
  // NOTE: we cannot use the Eigen:IncompleteLUT preconditioner, because its implementation is serial. Using it means
  // that BiCGSTAB will not run in parallel. It is faster without.

  solver.compute(A);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Eigen sparse solver failed at the compute step!");
  }

  vec_x = solver.solveWithGuess(b, guess);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Eigen sparse solver failed at the solve step!");
  }

  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error() << std::endl;

#pragma omp parallel for default(none) shared(arp, params, vec_x) collapse(2)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      // copy result into the wtd_T array:
      arp.wtd_T(x, y) = vec_x(arp.wtd_T.xyToI(x, y)) - arp.topo(x, y);
    }
  }
}

int update(Parameters& params, ArrayPack& arp, AppCtx& user_context, DMDA_Array_Pack& dmdapack) {
  PetscInt its;                // iterations for convergence
  SNESConvergedReason reason;  // Check convergence

  Eigen::initParallel();
  omp_set_num_threads(params.parallel_threads);
  Eigen::setNbThreads(params.parallel_threads);

  // compute any starting values needed for arrays
  set_starting_values(params, arp);

  // Get local array bounds
  const auto [xs, ys, xm, ym] = get_corners(user_context.da);

  // first_half(params, arp);
  //  values for storativity are reset each time; and recharge changes from one timestep to the next, so set these here
#pragma omp parallel for default(none) shared(arp, ys, ym, xs, xm, dmdapack, params) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.time_per_S[j][i] =
          params.deltat / updateEffectiveStorativity(arp.original_wtd(i, j), arp.wtd_T(i, j), arp.porosity(i, j));
      dmdapack.rech_vec[j][i] = add_recharge(arp.rech(i, j), arp.original_wtd(i, j), arp.porosity(i, j));
      dmdapack.head[j][i]     = arp.original_wtd(i, j) + arp.topo(i, j);
      dmdapack.guess[j][i]    = arp.wtd_T(i, j) + arp.topo(i, j);
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
  std::cout << "before solve 1 " << my_time() << std::endl;
  SNESSolve(user_context.snes, user_context.b, user_context.x);
  std::cout << "after solve 1 " << my_time() << std::endl;

  SNESGetIterationNumber(user_context.snes, &its);
  SNESGetConvergedReason(user_context.snes, &reason);

  PetscPrintf(
      PETSC_COMM_WORLD, "%s Number of nonlinear iterations = %" PetscInt_FMT "\n", SNESConvergedReasons[reason], its);

// recalculate the new storativity using the initial solve as an estimator for the final answer
#pragma omp parallel for default(none) shared(arp, ys, ym, xs, xm, dmdapack, params) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.time_per_S[j][i] =
          params.deltat / updateEffectiveStorativity(
                              arp.original_wtd(i, j), dmdapack.x[j][i] - dmdapack.topo_vec[j][i], arp.porosity(i, j));
    }
  }

// repeat the solve
#pragma omp parallel for default(none) shared(ys, ym, xs, xm, dmdapack) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      dmdapack.head[j][i]  = dmdapack.x[j][i];
      dmdapack.guess[j][i] = dmdapack.x[j][i];
    }
  }

  FormInitialGuess(&user_context, user_context.da, user_context.x);
  std::cout << "before solve 2 " << my_time() << std::endl;

  SNESSolve(user_context.snes, user_context.b, user_context.x);
  std::cout << "after solve 2 " << my_time() << std::endl;

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
  PetscScalar **x, **my_guess;

  DMDAVecGetArray(da, X, &x);
  DMDAVecGetArray(da, user_context->guess, &my_guess);

  const auto [xs, ys, xm, ym] = get_corners(da);

#pragma omp parallel for default(none) shared(my_guess, ys, ym, xs, xm, x) collapse(2)
  for (auto j = ys; j < ys + ym; j++) {
    for (auto i = xs; i < xs + xm; i++) {
      x[j][i] = my_guess[j][i];  // when land mask == 0, both topo and wtd have already been set to 0
                                 // elsewhere, so no need for another if statement here
    }
  }

  DMDAVecRestoreArray(da, X, &x);
  DMDAVecRestoreArray(da, user_context->guess, &my_guess);
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
  PetscScalar **time_over_storativity, **cellsize_ew_sq, **my_mask, **my_fdepth, **my_ksat, **my_topo, **my_rech,
      **my_T;

  /*
    Compute function over the locally owned part of the grid
 */
  PetscCall(DMDAVecGetArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecGetArray(da, user_context->time_per_S, &time_over_storativity));
  PetscCall(DMDAVecGetArray(da, user_context->cellsize_EW_squared, &cellsize_ew_sq));
  PetscCall(DMDAVecGetArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecGetArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecGetArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecGetArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecGetArray(da, user_context->T_vec, &my_T));

#pragma omp parallel for default(none) shared(info, my_T, x, my_topo, my_fdepth, my_ksat) collapse(2)
  for (auto j = info->ys; j < info->ys + info->ym; j++) {
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      my_T[j][i] = 1. / depthIntegratedTransmissivity(x[j][i] - my_topo[j][i], my_fdepth[j][i], my_ksat[j][i]);
    }
  }

#pragma omp parallel for default(none) \
    shared(info, cellsize_ew_sq, x, my_T, my_mask, my_rech, user_context, time_over_storativity, f)
  for (auto j = info->ys; j < info->ys + info->ym; j++) {
#pragma omp parallel for default(none) \
    shared(info, cellsize_ew_sq, my_mask, x, my_T, my_rech, user_context, j, time_over_storativity, f)
    for (auto i = info->xs; i < info->xs + info->xm; i++) {
      if (my_mask[j][i] == 0) {
        f[j][i] = 0;
      } else {
        const PetscScalar ux_E = (x[j][i + 1] - x[j][i]);
        const PetscScalar ux_W = (x[j][i] - x[j][i - 1]);
        const PetscScalar uy_N = (x[j + 1][i] - x[j][i]);
        const PetscScalar uy_S = (x[j][i] - x[j - 1][i]);
        const PetscScalar e_E  = 2. / (my_T[j][i] + my_T[j][i + 1]);
        const PetscScalar e_W  = 2. / (my_T[j][i] + my_T[j][i - 1]);
        const PetscScalar e_N  = 2. / (my_T[j][i] + my_T[j + 1][i]);
        const PetscScalar e_S  = 2. / (my_T[j][i] + my_T[j - 1][i]);

        const PetscScalar uxx = (e_W * ux_W - e_E * ux_E) / user_context->cellsize_NS_squared;
        const PetscScalar uyy = (e_S * uy_S - e_N * uy_N) / cellsize_ew_sq[j][i];

        f[j][i] = (uxx + uyy) * time_over_storativity[j][i] + x[j][i] -
                  my_rech[j][i];  // my_rech is converted to appropriate recharge for this timestep and starting water
                                  // table outside of the solve.
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(da, user_context->mask, &my_mask));
  PetscCall(DMDAVecRestoreArray(da, user_context->time_per_S, &time_over_storativity));
  PetscCall(DMDAVecRestoreArray(da, user_context->cellsize_EW_squared, &cellsize_ew_sq));
  PetscCall(DMDAVecRestoreArray(da, user_context->fdepth_vec, &my_fdepth));
  PetscCall(DMDAVecRestoreArray(da, user_context->ksat_vec, &my_ksat));
  PetscCall(DMDAVecRestoreArray(da, user_context->topo_vec, &my_topo));
  PetscCall(DMDAVecRestoreArray(da, user_context->rech_vec, &my_rech));
  PetscCall(DMDAVecRestoreArray(da, user_context->T_vec, &my_T));

  PetscLogFlops(info->xm * info->ym * (72.0));
  return 0;
}

}  // namespace FanDarcyGroundwater
