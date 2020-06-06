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
//
// This module is for packaging output quantities from RRTMGP based on spectral flux profiles
//    This implementation provides band-by-band flux profiles
//

#pragma once

#include "const_rrtmgpxx.h"
#include "mo_fluxes_byband_kernels.h"


class FluxesByband : public FluxesBroadband {
public:
  real3d bnd_flux_up;     // Band-by-band fluxes
  real3d bnd_flux_dn;     // (ncol, nlev, nband)
  real3d bnd_flux_net;    // Net (down - up)
  real3d bnd_flux_dn_dir; // Direct flux down


  template <class T>
  void reduce(real3d const &gpt_flux_up, real3d const &gpt_flux_dn, T const &spectral_disc, bool top_at_1, real3d const &gpt_flux_dn_dir=real3d()) {
    int ncol = size(gpt_flux_up, 1);
    int nlev = size(gpt_flux_up, 2);
    int ngpt = spectral_disc.get_ngpt();
    int nbnd = spectral_disc.get_nband();

    FluxesBroadband::reduce(gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir);

    if (size(gpt_flux_up,3) != ngpt) {
      stoprun("reduce: spectral discretization and g-point flux arrays have differing number of g-points");
    }
    // Check sizes of output arrays
    if (allocated(this->bnd_flux_up)) {
      if (size(this->band_flux_up,1) != ncol || size(this->band_flux_up,2) != nlev || size(this->band_flux_up,3) != nbnd) {
        stoprun("reduce: bnd_flux_up array incorrectly sized (can't compute net flux either)");
      }
    }
    if (allocated(this->bnd_flux_dn)) {
      if (size(this->band_flux_dn,1) != ncol || size(this->band_flux_dn,2) != nlev || size(this->band_flux_dn,3) != nbnd) {
        stoprun("reduce: bnd_flux_dn array incorrectly sized (can't compute net flux either)");
      }
    }
    if (allocated(this->bnd_flux_dir)) {
      if (size(this->band_flux_dir,1) != ncol || size(this->band_flux_dir,2) != nlev || size(this->band_flux_dir,3) != nbnd) {
        stoprun("reduce: bnd_flux_dir array incorrectly sized (can't compute net flux either)");
      }
    }
    if (allocated(this->bnd_flux_net)) {
      if (size(this->band_flux_net,1) != ncol || size(this->band_flux_net,2) != nlev || size(this->band_flux_net,3) != nbnd) {
        stoprun("reduce: bnd_flux_net array incorrectly sized (can't compute net flux either)");
      }
    }
    // Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if (allocated(this->bnd_flux_dn_dir) && ! allocated(gpt_flux_dn_dir)) {
      stoprun("reduce: requesting bnd_flux_dn_dir but direct flux hasn't been supplied");
    }

    int2d band_lims = spectral_disc.get_band_lims_gpoint();
    // Band-by-band fluxes
    // Up flux
    if (allocated(this->bnd_flux_up)) {
      sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_up,     this->bnd_flux_up    );
    }
    // Down flux
    if (allocated(this->bnd_flux_dn)) {
      sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn,     this->bnd_flux_dn    );
    }

    if (allocated(this->bnd_flux_dn_dir)) {
      sum_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn_dir, this->bnd_flux_dn_dir);
    }

    // Net flux
    if (allocated(this->bnd_flux_net)) {
      // Reuse down and up results if possible
      if (allocated(this->bnd_flux_dn) && allocated(this->bnd_flux_up)) {
        net_byband(ncol, nlev,       nbnd, this->bnd_flux_dn, this->bnd_flux_up, this->bnd_flux_net);
      } else {
        net_byband(ncol, nlev, ngpt, nbnd, band_lims, gpt_flux_dn, gpt_flux_up, this->bnd_flux_net);
      }
    }
  }
  
  
  // Are any fluxes desired from this set of g-point fluxes? We can tell because memory will be allocated for output
  bool are_desired() {
    return allocated(this->bnd_flux_up) || allocated(this->bnd_flux_dn) || allocated(this->bnd_flux_dn_dir) ||
           allocated(this->bnd_flux_net) || FluxesBroadband::are_desired();
  }


  void print_norms() const {
    FluxesBroadband::print_norms();
    if (allocated(bnd_flux_up    )) { std::cout << std::setprecision(16) << "bnd_flux_up    : " << sum(bnd_flux_up    ) << "\n"; }
    if (allocated(bnd_flux_dn    )) { std::cout << std::setprecision(16) << "bnd_flux_dn    : " << sum(bnd_flux_dn    ) << "\n"; }
    if (allocated(bnd_flux_net   )) { std::cout << std::setprecision(16) << "bnd_flux_net   : " << sum(bnd_flux_net   ) << "\n"; }
    if (allocated(bnd_flux_dn_dir)) { std::cout << std::setprecision(16) << "bnd_flux_dn_dir: " << sum(bnd_flux_dn_dir) << "\n"; }
  }
};

