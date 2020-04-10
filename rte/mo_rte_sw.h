// This code is part of Radiative Transfer for Energetics (RTE)
//
// Contacts: Robert Pincus and Eli Mlawer
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
// -------------------------------------------------------------------------------------------------
//
//  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
//    atmospheric optical properties on a spectral grid
//    information about vertical ordering
//    boundary conditions
//      solar zenith angle, spectrally-resolved incident colimated flux, surface albedos for direct and diffuse radiation
//    optionally, a boundary condition for incident diffuse radiation
//
// It is the user's responsibility to ensure that boundary conditions (incident fluxes, surface albedos) are on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------

namespace rte_sw {

  void rte_sw(OpticalProps2str const &atmos, bool top_at_1, real1d const &mu0, real2d const &inc_flux,
              real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, Fluxes &fluxes, real2d const &inc_flux_dif=real2d()) {
    real3d gpt_flux_up;
    real3d gpt_flux_dn;
    real3d gpt_flux_dir;
    real2d sfc_alb_dir_gpt;
    real2d sfc_alb_dif_gpt;
    int ncol  = atmos.get_ncol();
    int nlay  = atmos.get_nlay();
    int ngpt  = atmos.get_ngpt();
    int nband = atmos.get_nband();

    // Error checking -- consistency of sizes and validity of values
    if (! fluxes.are_desired()) { stoprun("rte_sw: no space allocated for fluxes"); }

    if (size(mu0,1) != ncol) { storun("rte_sw: mu0 inconsistently sized"); }
    if (anyLT(mu0,0._wp) || anyGT(mu0,1._wp)) { stoprun("rte_sw: one or more mu0 <= 0 or > 1"); }

    if (size(inc_flux,1) != ncol || size(inc_flux,2) != ngpt) { stoprun("rte_sw: inc_flux inconsistently sized"); }
    if (anyLT(inc_flux,0._wp)) { stoprun("rte_sw: one or more inc_flux < 0"); }
    if (allocated(inc_flux_dif)) {
      if (size(inc_flux_dif,1) != ncol || size(inc_flux_dif,2) != ngpt) { stoprun("rte_sw: inc_flux_dif inconsistently sized"); }
      if (anyLT(inc_flux_dif,0._wp)) { stoprun("rte_sw: one or more inc_flux_dif < 0"); }
    }

    if (size(sfc_alb_dir,1) != nband || size(sfc_alb_dir,2) != ncol) { stoprun("rte_sw: sfc_alb_dir inconsistently sized"); }
    if (anyLT(sfc_alb_dir,0._wp) || anyGT(sfc_alb_dir,1._wp)) { stoprun("rte_sw: sfc_alb_dir out of bounds [0,1]"); }
    if (size(sfc_alb_dif,1) != nband || size(sfc_alb_dif,2) != ncol) { stoprun("rte_sw: sfc_alb_dif inconsistently sized"); }
    if (anyLT(sfc_alb_dif,0._wp) || anyGT(sfc_alb_dif,1._wp)) { stoprun("rte_sw: sfc_alb_dif out of bounds [0,1]"); }

    gpt_flux_up  = real3d("gpt_flux_up" ,ncol, nlay+1, ngpt);
    gpt_flux_dn  = real3d("gpt_flux_dn" ,ncol, nlay+1, ngpt);
    gpt_flux_dir = real3d("gpt_flux_dir",ncol, nlay+1, ngpt);
    sfc_alb_dir_gpt = real2d("sfc_alb_dir_gpt",ncol, ngpt);
    sfc_alb_dif_gpt = real2d("sfc_alb_dif_gpt",ncol, ngpt);
    // Lower boundary condition -- expand surface albedos by band to gpoints
    //   and switch dimension ordering
    call expand_and_transpose(atmos, sfc_alb_dir, sfc_alb_dir_gpt);
    call expand_and_transpose(atmos, sfc_alb_dif, sfc_alb_dif_gpt);

    // Compute the radiative transfer...
    // Apply boundary conditions
    //   On input flux_dn is the diffuse component; the last action in each solver is to add
    //   direct and diffuse to represent the total, consistent with the LW
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, mu0, gpt_flux_dir);
    if (allocated(inc_flux_dif)) {
      call apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif,  gpt_flux_dn );
    } else {
      call apply_BC(ncol, nlay, ngpt, top_at_1,                gpt_flux_dn );
    }

    // two-stream calculation with scattering
    atmos.validate();
    sw_solver_2stream(ncol, nlay, ngpt, logical(top_at_1, wl), &
                      atmos.tau, atmos.ssa, atmos.g, mu0,      &
                      sfc_alb_dir_gpt, sfc_alb_dif_gpt,        &
                      gpt_flux_up, gpt_flux_dn, gpt_flux_dir);

    // ...and reduce spectral fluxes to desired output quantities
    fluxes.reduce(gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir);
  }


  // Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
  void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d &arr_out) {
    int ncol  = size(arr_in, 2)
    int nband = ops.get_nband()
    int ngpt  = ops.get_ngpt()
    int2d limits = ops.get_band_lims_gpoint();
    // for (int iband=1; iband <= nband; iband++) {
    //   for (int icol=1; icol <= ncol; icol++) {
    parallel_for( Bounds<2>(nband,ncol) , YAKL_LAMBDA (int iband, int icol) {
      for (int igpt=limits(1,iband); igpt <= limits(2,iband); igpt++) {
        arr_out(icol, igpt) = arr_in(iband,icol)
      }
    });
  }
}


