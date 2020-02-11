
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
YAKL_INLINE real interpolate3D(real1d const &scaling, real3d const &fmajor, real4d const &k, int igpt, int1d const &jeta,
                               int jtemp, int jpress, int ngpt, int neta, int npres, int ntemp) {
  // each code block is for a different reference temperature
  real ret = scaling(1) * ( fmajor(1,1,1) * k(igpt, jeta(1)  , jpress-1, jtemp  ) + 
                            fmajor(2,1,1) * k(igpt, jeta(1)+1, jpress-1, jtemp  ) + 
                            fmajor(1,2,1) * k(igpt, jeta(1)  , jpress  , jtemp  ) + 
                            fmajor(2,2,1) * k(igpt, jeta(1)+1, jpress  , jtemp  ) ) + 
             scaling(2) * ( fmajor(1,1,2) * k(igpt, jeta(2)  , jpress-1, jtemp+1) + 
                            fmajor(2,1,2) * k(igpt, jeta(2)+1, jpress-1, jtemp+1) + 
                            fmajor(1,2,2) * k(igpt, jeta(2)  , jpress  , jtemp+1) + 
                            fmajor(2,2,2) * k(igpt, jeta(2)+1, jpress  , jtemp+1) );
  return ret;
}



// ------------
//   This function returns a single value from a subset (in gpoint) of the k table
//
YAKL_INLINE real interpolate2D(real2d const &fminor, real3d const &k, int igpt, int1d const &jeta, int jtemp,
                               int ngpt, int neta, int ntemp) {
  return fminor(1,1) * k(igpt, jeta(1)  , jtemp  ) + 
         fminor(2,1) * k(igpt, jeta(1)+1, jtemp  ) + 
         fminor(1,2) * k(igpt, jeta(2)  , jtemp+1) + 
         fminor(2,2) * k(igpt, jeta(2)+1, jtemp+1);
}



// ----------------------------------------------------------
//
// One dimensional interpolation -- return all values along second table dimension
//
YAKL_INLINE void interpolate1D(real val, real offset, real delta, real2d const &table,
                               real2d &res, int tab_d1, int tab_d2) {
  real val0 = (val - offset) / delta;
  real frac = val0 - int(val0); // get fractional part
  int index = min(tab_d1-1, max(1, (int)(val0)+1)); // limit the index range
  for (int i=1; i<=tab_d2; i++) {
    res(i) = table(index,i) + frac * (table(index+1,i) - table(index,i));
  }
}



extern "C" void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta,
                                      int npres, int ntemp, int nPlanckTemp, real *tlay_p, real *tlev_p,
                                      real *tsfc_p, int sfc_lay, real *fmajor_p, int *jeta_p, bool *tropo_p,
                                      int *jtemp_p, int *jpress_p, int *gpoint_bands_p, int *band_lims_gpt_p,           
                                      real *pfracin_p, real temp_ref_min, real totplnk_delta, real *totplnk_p,
                                      int *gpoint_flavor_p, real *sfc_src_p, real *lay_src_p, real *lev_src_inc_p,
                                      real *lev_src_dec_p);



// ----------------------------------------------------------
//
// compute Rayleigh scattering optical depths
//
extern "C" void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta,
                                     int npres, int ntemp, int *gpoint_flavor_p, int *band_lims_gpt_p, 
                                     real *krayl_p, int idx_h2o, real *col_dry_p, real *col_gas_p,
                                     real *fminor_p, int *jeta_p, bool *tropo_p, int *jtemp_p,
                                     real *tau_rayleigh_p);



// ----------------------------------------------------------
//
// compute minor species optical depths
//
void gas_optical_depths_minor(int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2d &gpt_flv,
                              real3d &kminor, int2d &minor_limits_gpt, bool1d &minor_scales_with_density,
                              bool1d &scale_by_complement, int1d &idx_minor, int1d &idx_minor_scaling,
                              int1d &kminor_start, real2d &play, real2d &tlay, real3d &col_gas, real5d &fminor,
                              int4d &jeta, int2d &layer_limits, int2d &jtemp, real3d &tau);



// --------------------------------------------------------------------------------------
//
// compute minor species optical depths
//
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2d &gpoint_flavor, int2d &band_lims_gpt, real4d &kmajor,                         
                              real4d &col_mix, real6d &fmajor, int4d &jeta, bool2d &tropo, 
                              int2d &jtemp, int2d &jpress, real3d &tau);



// --------------------------------------------------------------------------------------
//
// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
//
extern "C" void compute_tau_absorption(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres,
                                       int ntemp, int nminorlower, int nminorklower, int nminorupper, int nminorkupper,          
                                       int idx_h2o, int *gpoint_flavor_p, int *band_lims_gpt_p, real *kmajor_p,
                                       real *kminor_lower_p, real *kminor_upper_p, int *minor_limits_gpt_lower_p,
                                       int *minor_limits_gpt_upper_p, bool *minor_scales_with_density_lower_p,    
                                       bool *minor_scales_with_density_upper_p, bool *scale_by_complement_lower_p,          
                                       bool *scale_by_complement_upper_p, int *idx_minor_lower_p,                    
                                       int *idx_minor_upper_p, int *idx_minor_scaling_lower_p,            
                                       int *idx_minor_scaling_upper_p, int *kminor_start_lower_p,                 
                                       int *kminor_start_upper_p, bool *tropo_p, real *col_mix_p, real *fmajor_p,
                                       real *fminor_p, real *play_p, real *tlay_p, real *col_gas_p, int *jeta_p,
                                       int *jtemp_p, int *jpress_p, real *tau_p);


// Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
//   using Rayleigh scattering phase function
extern "C" void combine_and_reorder_nstr(int ncol, int nlay, int ngpt, int nmom, real *tau_abs_p, real *tau_rayleigh_p, real *tau_p, real *ssa_p, real *p_p);




