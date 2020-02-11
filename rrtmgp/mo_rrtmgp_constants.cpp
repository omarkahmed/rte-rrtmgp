
#include "mo_rrtmgp_constants.h"


extern "C" void init_constants(real gravity, real mol_weight_dry_air, real heat_capacity_dry_air) {
  grav   = gravity;
  m_dry  = mol_weight_dry_air;
  cp_dry = heat_capacity_dry_air;
}



