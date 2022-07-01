#include "parameters.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

struct AppCtx {
  PetscReal timestep;
  PetscReal cellsize_NS;
  // SNES snes       = nullptr;
  DM da           = nullptr;
  Vec x           = nullptr;  // Solution vector
  Vec xdot        = nullptr;  // Solution vector
  Vec b           = nullptr;  // RHS vector
  Vec S           = nullptr;
  Vec cellsize_EW = nullptr;
  Vec fdepth_vec  = nullptr;
  Vec ksat_vec    = nullptr;
  Vec mask        = nullptr;
  Vec porosity    = nullptr;
  Vec h           = nullptr;
  Vec topo_vec    = nullptr;
  Vec rech_vec    = nullptr;
  Vec cell_area   = nullptr;
  PetscReal lidvelocity, prandtl, grashof; /* physical parameters */
  PetscBool parabolic;   /* allow a transient term corresponding roughly to artificial compressibility */
  PetscReal cfl_initial; /* CFL for first time step */

  ~AppCtx() {
    //  SNESDestroy(&snes);
    DMDestroy(&da);
    VecDestroy(&x);
    VecDestroy(&xdot);
    VecDestroy(&b);
    VecDestroy(&S);
    VecDestroy(&cellsize_EW);
    VecDestroy(&fdepth_vec);
    VecDestroy(&ksat_vec);
    VecDestroy(&mask);
    VecDestroy(&porosity);
    VecDestroy(&h);
    VecDestroy(&topo_vec);
    VecDestroy(&rech_vec);
    VecDestroy(&cell_area);
  }

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &b);
    VecDuplicate(x, &xdot);
    VecDuplicate(x, &S);
    VecDuplicate(x, &cellsize_EW);
    VecDuplicate(x, &fdepth_vec);
    VecDuplicate(x, &ksat_vec);
    VecDuplicate(x, &mask);
    VecDuplicate(x, &porosity);
    VecDuplicate(x, &h);
    VecDuplicate(x, &topo_vec);
    VecDuplicate(x, &rech_vec);
    VecDuplicate(x, &cell_area);
  }
};

void InitialiseSNES(AppCtx& user_context, Parameters& params);
