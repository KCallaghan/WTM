#include <experimental/source_location>

struct DMDA_Array_Pack {
  PetscScalar** x                   = nullptr;
  PetscScalar** cellsize_EW_squared = nullptr;
  PetscScalar** fdepth_vec          = nullptr;
  PetscScalar** ksat_vec            = nullptr;
  PetscScalar** mask                = nullptr;
  PetscScalar** topo_vec            = nullptr;
  PetscScalar** rech_vec            = nullptr;
  PetscScalar** T_vec               = nullptr;
  PetscScalar** head                = nullptr;
  PetscScalar** porosity_vec        = nullptr;
  PetscScalar** starting_wtd        = nullptr;
  const AppCtx* context             = nullptr;

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    DMDAVecGetArray(user.da, user.x, &x);
    DMDAVecGetArray(user.da, user.cellsize_EW_squared, &cellsize_EW_squared);
    DMDAVecGetArray(user.da, user.fdepth_vec, &fdepth_vec);
    DMDAVecGetArray(user.da, user.ksat_vec, &ksat_vec);
    DMDAVecGetArray(user.da, user.mask, &mask);
    DMDAVecGetArray(user.da, user.topo_vec, &topo_vec);
    DMDAVecGetArray(user.da, user.rech_vec, &rech_vec);
    DMDAVecGetArray(user.da, user.T_vec, &T_vec);
    DMDAVecGetArray(user.da, user.head, &head);
    DMDAVecGetArray(user.da, user.porosity_vec, &porosity_vec);
    DMDAVecGetArray(user.da, user.starting_wtd, &starting_wtd);
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    DMDAVecRestoreArray(context->da, context->x, &x);
    DMDAVecRestoreArray(context->da, context->cellsize_EW_squared, &cellsize_EW_squared);
    DMDAVecRestoreArray(context->da, context->fdepth_vec, &fdepth_vec);
    DMDAVecRestoreArray(context->da, context->ksat_vec, &ksat_vec);
    DMDAVecRestoreArray(context->da, context->mask, &mask);
    DMDAVecRestoreArray(context->da, context->topo_vec, &topo_vec);
    DMDAVecRestoreArray(context->da, context->rech_vec, &rech_vec);
    DMDAVecRestoreArray(context->da, context->T_vec, &T_vec);
    DMDAVecRestoreArray(context->da, context->head, &head);
    DMDAVecRestoreArray(context->da, context->porosity_vec, &porosity_vec);
    DMDAVecRestoreArray(context->da, context->starting_wtd, &starting_wtd);
    context = nullptr;
  }
};

void populate_DMDA_array_pack(AppCtx& user_context, ArrayPack& arp, DMDA_Array_Pack& dmdapack);
