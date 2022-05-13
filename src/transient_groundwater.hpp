#pragma once

#include "ArrayPack.hpp"
#include "parameters.hpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscerror.h>
#include <petscsnes.h>

#include <richdem/common/Array2D.hpp>

namespace FanDarcyGroundwater {

void update(Parameters& params, ArrayPack& arp);

}
