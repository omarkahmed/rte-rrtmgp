
#pragma once

#include "const.h"


// Spectral reduction over all points
extern "C" void sum_broadband(int ncol, int nlev, int ngpt, real *spectral_flux_p, real *broadband_flux_p);


// Net flux: Spectral reduction over all points
extern "C" void net_broadband_full(int ncol, int nlev, int ngpt, real *spectral_flux_dn_p, real *spectral_flux_up_p, real *broadband_flux_net_p);


// Net flux when bradband flux up and down are already available
extern "C" void net_broadband_precalc(int ncol, int nlev, real *flux_dn_p, real *flux_up_p, real *broadband_flux_net_p);


