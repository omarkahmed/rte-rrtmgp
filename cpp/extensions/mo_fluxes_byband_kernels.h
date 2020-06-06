// This code is part of Radiative Transfer for Energetics (RTE)
//
// Eli Mlawer and Robert Pincus
// Andre Wehe and Jennifer Delamere
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
//
// Kernels for computing byband fluxes by summing over all elements in the spectral dimension
//

// Spectral reduction over all points
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &band_lims, real3d const &spectral_flux, real3d &byband_flux) {
  // do ibnd = 1, nbnd
  //   do ilev = 1, nlev
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA(int ibnd, int ilev, int icol) {
    byband_flux(icol, ilev, ibnd) =  spectral_flux(icol, ilev, band_lims(1, ibnd));
    for (int igpt = band_lims(1,ibnd)+1 , igpt <= band_lims(2,ibnd) , igpt++) {
      byband_flux(icol, ilev, ibnd) = byband_flux(icol, ilev, ibnd) + spectral_flux(icol, ilev, igpt);
    }
  });
}


// Net flux: Spectral reduction over all points
void net_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &band_lims, real3d const &spectral_flux_dn,
                real3d const &spectral_flux_up, real3d &byband_flux_net) {
  // do ibnd = 1, nbnd
  //   do ilev = 1, nlev
  //     do icol = 1, ncol
  parallel_for( Bounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA (int ibnd, int ilev, int icol) {
    int igpt = band_lims(1,ibnd);
    byband_flux_net(icol, ilev, ibnd) = spectral_flux_dn(icol, ilev, igpt) - spectral_flux_up(icol, ilev, igpt);
    for (igpt = band_lims(1,ibnd)+1 ; igpt <= band_lims(2,ibnd) ; igpt++) {
      byband_flux_net(icol, ilev, ibnd) = byband_flux_net (icol, ilev, ibnd) +
                                          spectral_flux_dn(icol, ilev, igpt) -
                                          spectral_flux_up(icol, ilev, igpt);
    }
  });
}


void net_byband(int ncol, int nlev, int nbnd, real3d const &byband_flux_dn, real3d const &byband_flux_up, real3d &byband_flux_net) {
  parallel_for( Bounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA (int ibnd, int ilev, int icol) {
    byband_flux_net(icol,ilev,ibnd) = byband_flux_dn(icol,ilev,ibnd) - byband_flux_up(icol,ilev,ibnd)
  });
}

