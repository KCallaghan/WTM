//#include "CreateSNES.hpp"

void InitialiseSNES(AppCtx& user_context, Parameters& params) {
  SNESCreate(PETSC_COMM_WORLD, &user_context.snes);

  user_context.cellsize_NS = params.cellsize_n_s_metres;
  user_context.timestep    = params.deltat;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMDACreate2d(
      PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE,
      DM_BOUNDARY_NONE,
      DMDA_STENCIL_STAR,
      params.ncells_x,
      params.ncells_y,
      PETSC_DECIDE,
      PETSC_DECIDE,
      1,
      1,
      nullptr,
      nullptr,
      &user_context.da);
  DMSetFromOptions(user_context.da);
  DMSetUp(user_context.da);
}
