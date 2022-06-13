#include <experimental/source_location>

// void PETSC_CHECK(
//     const PetscErrorCode err,
//     const std::experimental::source_location location = std::experimental::source_location::current()) {
//   if (err) {
//     throw std::runtime_error(
//         "Petsc exception: " + std::to_string(err) + " at " + location.file_name() + ":" +
//         std::to_string(location.line()));
//   }
// }

struct DMDA_Array_Pack {
  PetscScalar** x           = nullptr;
  PetscScalar** S           = nullptr;
  PetscScalar** cellsize_EW = nullptr;
  PetscScalar** fdepth_vec  = nullptr;
  PetscScalar** ksat_vec    = nullptr;
  PetscScalar** mask        = nullptr;
  PetscScalar** porosity    = nullptr;
  PetscScalar** h           = nullptr;
  PetscScalar** topo_vec    = nullptr;
  PetscScalar** rech_vec    = nullptr;
  const AppCtx* context     = nullptr;

  DMDA_Array_Pack(const AppCtx& user) {
    assert(!context);  // Make sure we're not already initialized
    context = &user;
    DMDAVecGetArray(user.da, user.x, &x);
    DMDAVecGetArray(user.da, user.S, &S);
    DMDAVecGetArray(user.da, user.cellsize_EW, &cellsize_EW);
    DMDAVecGetArray(user.da, user.fdepth_vec, &fdepth_vec);
    DMDAVecGetArray(user.da, user.ksat_vec, &ksat_vec);
    DMDAVecGetArray(user.da, user.mask, &mask);
    DMDAVecGetArray(user.da, user.porosity, &porosity);
    DMDAVecGetArray(user.da, user.h, &h);
    DMDAVecGetArray(user.da, user.topo_vec, &topo_vec);
    DMDAVecGetArray(user.da, user.rech_vec, &rech_vec);
  }

  void release() {
    assert(context);  // Make sure we are already initialized
    DMDAVecRestoreArray(context->da, context->x, &x);
    DMDAVecRestoreArray(context->da, context->S, &S);
    DMDAVecRestoreArray(context->da, context->cellsize_EW, &cellsize_EW);
    DMDAVecRestoreArray(context->da, context->fdepth_vec, &fdepth_vec);
    DMDAVecRestoreArray(context->da, context->ksat_vec, &ksat_vec);
    DMDAVecRestoreArray(context->da, context->mask, &mask);
    DMDAVecRestoreArray(context->da, context->porosity, &porosity);
    DMDAVecRestoreArray(context->da, context->h, &h);
    DMDAVecRestoreArray(context->da, context->topo_vec, &topo_vec);
    DMDAVecRestoreArray(context->da, context->rech_vec, &rech_vec);
    context = nullptr;
  }
};

void populate_DMDA_array_pack(AppCtx& user_context, ArrayPack& arp, DMDA_Array_Pack& dmdapack);
