
#pragma once
#include "const.h"
#include "YAKL.h"
#include "YAKL_netcdf.h"
#include "mo_gas_concentrations.h"

void read_atmos(std::string input_file, realHost2d &p_lay, realHost2d &t_lay, realHost2d &p_lev, realHost2d &t_lev,
                GasConcs &gas_concs, realHost2d &col_dry, int ncol);
