
#pragma once

#include "const.h"

// This code is part of
// RRTM for GCM Applications - Parallel (RRTMGP)
//
// Eli Mlawer and Robert Pincus
// Andre Wehe and Jennifer Delamere
// email:  rrtmgp@aer.com
//
// Copyright 2015,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
//
// Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths,
//   source functions.




// Compute interpolation coefficients
// for calculations of major optical depths, minor optical depths, Rayleigh,
// and Planck fractions
extern "C" void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp,
                              int *flavor_p, real *press_ref_log_p, real *temp_ref_p, real press_ref_log_delta,
                              real temp_ref_min, real temp_ref_delta, real press_ref_trop_log,
                              real *vmr_ref_p, real *play_p, real *tlay_p, real *col_gas_p, int *jtemp_p,
                              real *fmajor_p, real *fminor_p, real *col_mix_p, bool *tropo_p, int *jeta_p,
                              int *jpress_p);


// ----------------------------------------------------------
//
// Combine absoprtion and Rayleigh optical depths for total tau, ssa, g
//
extern "C" void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real *tau_abs_p, real *tau_rayleigh_p,
                                         real *tau_p,  real *ssa_p, real *g_p);



//  ----------------------------------------------------------
//  interpolation in temperature, pressure, and eta
extern "C" real interpolate3D(real *scaling_p, real *fmajor_p, real *k_p, int igpt, int *jeta_p,
                              int jtemp, int jpress, int ngpt, int neta, int npres, int ntemp);
real interpolate3D(real1d scaling, real3d fmajor, real4d k, int igpt, int1d jeta,
                   int jtemp, int jpress, int ngpt, int neta, int npres, int ntemp);



// ------------
//   This function returns a single value from a subset (in gpoint) of the k table
//
extern "C" real interpolate2D(real *fminor_p, real *k_p, int igpt, int *jeta_p, int jtemp,
                              int ngpt, int neta, int ntemp);
real interpolate2D(real2d fminor, real3d k, int igpt, int1d jeta, int jtemp,
                   int ngpt, int neta, int ntemp);



// ----------------------------------------------------------
//
// One dimensional interpolation -- return all values along second table dimension
//
extern "C" void interpolate1D(real val, real offset, real delta, real *table_p,
                              real *res_p, int tab_d1, int tab_d2);
void interpolate1D(real val, real offset, real delta, real2d table,
                   real2d res, int tab_d1, int tab_d2);




extern "C" void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta,
                                      int npres, int ntemp, int nPlanckTemp, real *tlay_p, real *tlev_p,
                                      real *tsfc_p, int sfc_lay, real *fmajor_p, int *jeta_p, bool *tropo_p,
                                      int *jtemp_p, int *jpress_p, int *gpoint_bands_p, int *band_lims_gpt_p,           
                                      real *pfracin_p, real temp_ref_min, real totplnk_delta, real *totplnk_p,
                                      int *gpoint_flavor_p, real *sfc_src_p, real *lay_src_p, real *lev_src_inc_p,
                                      real *lev_src_dec_p);
