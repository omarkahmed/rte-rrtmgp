
#pragma once

#include "const_rrtmgpxx.h"
#include "mo_gas_concentrations.h"
#include "mo_gas_optics_rrtmgp.h"
#include <string>

void load_and_init(GasOpticsRRTMGP &kdist, std::string filename, GasConcs const &available_gases);

