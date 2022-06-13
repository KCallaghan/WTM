#pragma once

#include "ArrayPack.hpp"
#include "CreateSNES.hpp"
//#include "DMDA_array_pack.hpp"

#include <richdem/common/Array2D.hpp>

namespace FanDarcyGroundwater {

int update(Parameters& params, ArrayPack& arp, AppCtx& user_context);

}
