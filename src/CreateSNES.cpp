//#include "CreateSNES.hpp"

void InitialiseSNES(AppCtx& user_context, Parameters& params) {
  // SNESCreate(PETSC_COMM_WORLD, &user_context.snes);
  TSCreate(PETSC_COMM_WORLD, &user_context.ts);
  TSSetType(user_context.ts, TSBDF);
  TSBDFSetOrder(user_context.ts, 1);  // BDF1 = Backward Euler

  user_context.cellsize_NS = params.cellsize_n_s_metres;
  user_context.maxtime     = params.deltat;
  user_context.timestep    = params.petsc_timestep;

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

  user_context.make_global_vectors();

  DMSetApplicationContext(user_context.da, &user_context);
  //  SNESSetDM(user_context.snes, user_context.da);
  TSSetDM(user_context.ts, user_context.da);

  //  SNESSetFromOptions(user_context.snes);

  //    PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"TS conservation example","");
  // PetscOptionsEnum("-var","Variable formulation",NULL,VarModes,(PetscEnum)var,(PetscEnum*)&var,NULL);
  // PetscOptionsEnd();
}
