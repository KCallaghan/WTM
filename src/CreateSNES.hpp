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
  Vec cellsize_EW_squared = nullptr;
  Vec fdepth_vec          = nullptr;
  Vec ksat_vec            = nullptr;
  Vec mask                = nullptr;
  Vec topo_vec            = nullptr;
  Vec rech_vec            = nullptr;
  Vec porosity_vec        = nullptr;
  Vec starting_wtd        = nullptr;

  // Local ghost vectors for fields accessed at neighbor indices in FormFunctionLocal
  Vec topo_local   = nullptr;
  Vec fdepth_local = nullptr;
  Vec ksat_local   = nullptr;
  Vec T_local      = nullptr;  // scratch: 1/T, computed over ghost range each F eval

  // Extract global vectors from DM; then duplicate for remaining
  // vectors that are the same types
  void make_global_vectors() {
    DMCreateGlobalVector(da, &x);
    VecDuplicate(x, &b);
    VecDuplicate(x, &cellsize_EW_squared);
    VecDuplicate(x, &fdepth_vec);
    VecDuplicate(x, &ksat_vec);
    VecDuplicate(x, &mask);
    VecDuplicate(x, &topo_vec);
    VecDuplicate(x, &rech_vec);
    VecDuplicate(x, &porosity_vec);
    VecDuplicate(x, &starting_wtd);
  }

  void make_local_vectors() {
    DMCreateLocalVector(da, &topo_local);
    DMCreateLocalVector(da, &fdepth_local);
    DMCreateLocalVector(da, &ksat_local);
    DMCreateLocalVector(da, &T_local);
  }
};

void InitialiseSNES(AppCtx& user_context, Parameters& params);
