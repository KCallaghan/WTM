#include "parameters.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

struct AppCtx {
  PetscReal cellsize_NS_squared;
  SNES snes               = nullptr;
  DM da                   = nullptr;
  Vec x                   = nullptr;  // Solution vector
  Vec b                   = nullptr;  // RHS vector
  Vec time_per_S          = nullptr;
  Vec cellsize_EW_squared = nullptr;
  Vec fdepth_vec          = nullptr;
  Vec ksat_vec            = nullptr;
  Vec mask                = nullptr;
  Vec topo_vec            = nullptr;
  Vec rech_vec            = nullptr;
  Vec T_vec               = nullptr;

  ~AppCtx() {
    SNESDestroy(&snes);
    DMDestroy(&da);
    VecDestroy(&x);
    VecDestroy(&b);
    VecDestroy(&time_per_S);
    VecDestroy(&cellsize_EW_squared);
    VecDestroy(&fdepth_vec);
    VecDestroy(&ksat_vec);
    VecDestroy(&mask);
    VecDestroy(&topo_vec);
    VecDestroy(&rech_vec);
    VecDestroy(&T_vec);
  }

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &b);
    VecDuplicate(x, &time_per_S);
    VecDuplicate(x, &cellsize_EW_squared);
    VecDuplicate(x, &fdepth_vec);
    VecDuplicate(x, &ksat_vec);
    VecDuplicate(x, &mask);
    VecDuplicate(x, &topo_vec);
    VecDuplicate(x, &rech_vec);
    VecDuplicate(x, &T_vec);
  }
};

void InitialiseSNES(AppCtx& user_context, Parameters& params);
