#include "transient_groundwater.hpp"
#include "add_recharge.hpp"
#include "update_effective_storativity.hpp"

#include <omp.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>  //obtained on Linux using apt install libeigen3-dev. Make sure this points to the right place to include.

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;  // declares a row-major sparse matrix type of double

constexpr double seconds_in_a_year = 31536000.;

///////////////////////
// PRIVATE FUNCTIONS //
///////////////////////

namespace FanDarcyGroundwater {

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

  std::cout << "do b and the guess" << std::endl;
#pragma omp parallel for default(none) shared(arp, b, guess, params) collapse(2)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      b(arp.wtd_T.xyToI(x, y)) =
          arp.wtd(x, y) + arp.topo(x, y);  // wtd is 0 in ocean cells and topo is 0 in ocean cells, so no need to
                                           // differentiate between ocean vs land.
      guess(arp.wtd_T.xyToI(x, y)) = arp.wtd_T(x, y) + arp.topo(x, y);
    }
  }

  std::cout << "do the matrix" << std::endl;
#pragma omp parallel for collapse(2) default(none) shared(arp, params, A)
  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      // The row and column that the current cell will be stored in in matrix A.
      // This should go up monotonically, i.e. [0,0]; [1,1]; [2,2]; etc.
      // All of the N,E,S,W directions should be in the same row, but the column will differ.
      const auto main_loc = arp.wtd_T.xyToI(x, y);

      int d_w = (x == 0) ? 0 : -1;
      int d_e = (x == params.ncells_x - 1) ? 0 : 1;
      int d_s = (y == 0) ? 0 : -1;
      int d_n = (y == params.ncells_y - 1) ? 0 : 1;

      double T_N =
          2. /
          (1. / (depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y))) +
           1. / (depthIntegratedTransmissivity(arp.wtd(x, y + d_n), arp.fdepth(x, y + d_n), arp.ksat(x, y + d_n))));
      double T_S =
          2. /
          (1. / (depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y))) +
           1. / (depthIntegratedTransmissivity(arp.wtd(x, y + d_s), arp.fdepth(x, y + d_s), arp.ksat(x, y + d_s))));
      double T_E =
          2. /
          (1. / (depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y))) +
           1. / (depthIntegratedTransmissivity(arp.wtd(x + d_e, y), arp.fdepth(x + d_e, y), arp.ksat(x + d_e, y))));
      double T_W =
          2. /
          (1. / (depthIntegratedTransmissivity(arp.wtd(x, y), arp.fdepth(x, y), arp.ksat(x, y))) +
           1. / (depthIntegratedTransmissivity(arp.wtd(x + d_w, y), arp.fdepth(x + d_w, y), arp.ksat(x + d_w, y))));

      if (x != 0) {
        // Do the North diagonal. Offset by -1. When x == 0, there is no north diagonal.
        A.insert(main_loc, main_loc - 1) =
            (-T_N * params.deltat) /
            (arp.effective_storativity(x, y) * params.cellsize_n_s_metres * params.cellsize_n_s_metres);
      }

      if (y != 0) {
        // Next is the West diagonal. Opposite of the East. Located at (i,j-params.ncells_x). When y == 0, there is no
        // west diagonal.
        A.insert(main_loc, main_loc - params.ncells_x) =
            (-T_W * params.deltat) /
            (arp.effective_storativity(x, y) * arp.cellsize_e_w_metres[y] * arp.cellsize_e_w_metres[y]);
      }

      // major diagonal:
      A.insert(main_loc, main_loc) =
          ((((T_E + T_W) * params.deltat) /
            (arp.effective_storativity(x, y) * arp.cellsize_e_w_metres[y] * arp.cellsize_e_w_metres[y])) +
           (((T_N + T_S) * params.deltat) /
            (arp.effective_storativity(x, y) * params.cellsize_n_s_metres * params.cellsize_n_s_metres)) +
           1);

      if (y != params.ncells_y - 1) {
        // Now do the East diagonal. The East location is at (i,j+params.ncells_x). When y == params.ncells_y -1,
        // there is no east diagonal.
        A.insert(main_loc, main_loc + params.ncells_x) =
            (-T_E * params.deltat) /
            (arp.effective_storativity(x, y) * arp.cellsize_e_w_metres[y] * arp.cellsize_e_w_metres[y]);
      }

      if (x != params.ncells_x - 1) {
        // Do the South diagonal, offset by +1. When x == params.ncells_x, there is no south diagonal.
        A.insert(main_loc, main_loc + 1) =
            (-T_S * params.deltat) /
            (arp.effective_storativity(x, y) * params.cellsize_n_s_metres * params.cellsize_n_s_metres);
      }
    }
  }
  std::cout << "makeCompressed" << std::endl;
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
      // copy result into the wtd array:
      arp.wtd(x, y) = vec_x(arp.wtd_T.xyToI(x, y)) - arp.topo(x, y);
    }
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
        arp.wtd_T(x, y)                 = 0.;
        arp.original_wtd(x, y)          = 0.;
        arp.effective_storativity(x, y) = 1.;
      } else {
        arp.original_wtd(x, y) = arp.wtd(x, y);
        // Apply the first half of the recharge to the water-table depth grid (wtd)
        // Its clone (wtd_T) is used and updated in the Picard iteration
        // use regular porosity for adding recharge since this checks
        // for underground space within add_recharge.
        arp.wtd(x, y) += add_recharge(arp.rech(x, y), arp.wtd(x, y), arp.porosity(x, y), arp.cell_area[y], 1, arp);

        arp.wtd_T(x, y) = arp.wtd(x, y);

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

void update(Parameters& params, ArrayPack& arp) {
  Eigen::initParallel();
  omp_set_num_threads(params.parallel_threads);
  Eigen::setNbThreads(params.parallel_threads);

  // set starting values and arrays for the groundwater calculation
  set_starting_values(params, arp);

  // Picard iteration through solver:
  // Repeat the first solve params.picard_iterations number of times.
  // This helps us to get a more accurate halfway transmissivity
  // to use in the second_half.
  for (int continue_picard = 0; continue_picard < params.picard_iterations; continue_picard++) {
    // Run the first half solver: solve for water table at half of the total deltat.
    std::cout << "p first_half" << std::endl;
    first_half(params, arp);

// update the effective storativity to use during the next iteration of the first half:
#pragma omp parallel for default(none) shared(arp, params) collapse(2)
    for (int y = 0; y < params.ncells_y; y++) {
      for (int x = 0; x < params.ncells_x; x++) {
        if (arp.land_mask(x, y) != 0.f) {
          arp.effective_storativity(x, y) = updateEffectiveStorativity(
              arp.original_wtd(x, y), arp.wtd_T(x, y), arp.porosity(x, y), arp.effective_storativity(x, y));
        }
      }
    }
  }

  for (int y = 0; y < params.ncells_y; y++) {
    for (int x = 0; x < params.ncells_x; x++) {
      if (arp.land_mask(x, y) != 0) {
        continue;
      }
      // count up any water lost to the ocean so that we can compute a water balance
      if (arp.wtd(x, y) > 0) {
        arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y];
      } else {
        arp.total_loss_to_ocean += arp.wtd(x, y) * arp.cell_area[y] * arp.porosity(x, y);
      }
      arp.wtd(x, y) = 0.;  // reset ocean water tables to 0.
    }
  }
}

}  // namespace FanDarcyGroundwater
