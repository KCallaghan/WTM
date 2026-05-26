#include <experimental/source_location>

struct DMDA_Array_Pack {
  PetscScalar** x                   = nullptr;
  PetscScalar** cellsize_EW_squared = nullptr;
  PetscScalar** mask                = nullptr;
  PetscScalar** rech_vec            = nullptr;
  PetscScalar** porosity_vec        = nullptr;
  PetscScalar** starting_wtd        = nullptr;
  const AppCtx* context             = nullptr;

  // topo_vec, fdepth_vec, ksat_vec are intentionally NOT held here.
  // They are scattered to AppCtx local ghost vectors before each solve so that
  // FormFunctionLocal can safely access neighbor indices across MPI boundaries.

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    DMDAVecGetArray(user.da, user.x, &x);
    DMDAVecGetArray(user.da, user.cellsize_EW_squared, &cellsize_EW_squared);
    DMDAVecGetArray(user.da, user.mask, &mask);
    DMDAVecGetArray(user.da, user.rech_vec, &rech_vec);
    DMDAVecGetArray(user.da, user.porosity_vec, &porosity_vec);
    DMDAVecGetArray(user.da, user.starting_wtd, &starting_wtd);
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    DMDAVecRestoreArray(context->da, context->x, &x);
    DMDAVecRestoreArray(context->da, context->cellsize_EW_squared, &cellsize_EW_squared);
    DMDAVecRestoreArray(context->da, context->mask, &mask);
    DMDAVecRestoreArray(context->da, context->rech_vec, &rech_vec);
    DMDAVecRestoreArray(context->da, context->porosity_vec, &porosity_vec);
    DMDAVecRestoreArray(context->da, context->starting_wtd, &starting_wtd);
    context = nullptr;
  }
};

void populate_DMDA_array_pack(AppCtx& user_context, ArrayPack& arp, DMDA_Array_Pack& dmdapack);
void scatter_static_fields(AppCtx& user_context, ArrayPack& arp);
