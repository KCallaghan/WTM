#include "parameters.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

struct AppCtx {
  PetscReal cellsize_NS_squared;
  PetscReal deltat;
  SNES snes               = nullptr;
  DM da                   = nullptr;
  Vec x                   = nullptr;  // Solution vector
  Vec b                   = nullptr;  // RHS vector
  Vec storativity_vec     = nullptr;
  Vec cellsize_EW_squared = nullptr;
  Vec fdepth_vec          = nullptr;
  Vec ksat_vec            = nullptr;
  Vec mask                = nullptr;
  Vec topo_vec            = nullptr;
  Vec rech_vec            = nullptr;
  Vec T_vec               = nullptr;
  Vec head                = nullptr;
  Vec porosity_vec        = nullptr;

  ~AppCtx() {
    SNESDestroy(&snes);
    DMDestroy(&da);
    VecDestroy(&x);
    VecDestroy(&b);
    VecDestroy(&storativity_vec);
    VecDestroy(&cellsize_EW_squared);
    VecDestroy(&fdepth_vec);
    VecDestroy(&ksat_vec);
    VecDestroy(&mask);
    VecDestroy(&topo_vec);
    VecDestroy(&rech_vec);
    VecDestroy(&T_vec);
    VecDestroy(&head);
    VecDestroy(&porosity_vec);
  }

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &b);
    VecDuplicate(x, &storativity_vec);
    VecDuplicate(x, &cellsize_EW_squared);
    VecDuplicate(x, &fdepth_vec);
    VecDuplicate(x, &ksat_vec);
    VecDuplicate(x, &mask);
    VecDuplicate(x, &topo_vec);
    VecDuplicate(x, &rech_vec);
    VecDuplicate(x, &T_vec);
    VecDuplicate(x, &head);
    VecDuplicate(x, &porosity_vec);
  }
};

void InitialiseSNES(AppCtx& user_context, Parameters& params);
