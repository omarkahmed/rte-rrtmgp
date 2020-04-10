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
//    atmospheric optical properties, spectrally-resolved
//    information about vertical ordering
//    internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
//    boundary conditions: surface emissivity defined per band
//    optionally, a boundary condition for incident diffuse radiation
//    optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected
//
// If optical properties are supplied via class ty_optical_props_1scl (absorption optical thickenss only)
//    then an emission/absorption solver is called
//    If optical properties are supplied via class ty_optical_props_2str fluxes are computed via
//    two-stream calculations and adding.
//
// It is the user's responsibility to ensure that emissivity is on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------

namespace rte_lw {

  // Interface using only optical properties and source functions as inputs; fluxes as outputs.
  void rte_lw(OpticalProps1scl const &optical_props, bool top_at_1, SourceFuncLW const &sources, real2d const &sfc_emis,
              Fluxes &fluxes, real2d const &inc_flux=real2d(), int n_gauss_angles=-1) {
    real3d gpt_flux_up;
    real3d gpt_flux_dn;
    real2d sfc_emis_gpt;

    // Weights and angle secants for first order (k=1) Gaussian quadrature.
    //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
    //   after Abramowitz & Stegun 1972, page 921
    int constexpr max_gauss_pts = 4;
    real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
    gauss_Ds(1,1) = 1.66_wp      ; gauss_Ds(2,1) =         0._wp; gauss_Ds(3,1) =         0._wp; gauss_Ds(4,1) =         0._wp;
    gauss_Ds(1,2) = 1.18350343_wp; gauss_Ds(2,2) = 2.81649655_wp; gauss_Ds(3,2) =         0._wp; gauss_Ds(4,2) =         0._wp;
    gauss_Ds(1,3) = 1.09719858_wp; gauss_Ds(2,3) = 1.69338507_wp; gauss_Ds(3,3) = 4.70941630_wp; gauss_Ds(4,3) =         0._wp;
    gauss_Ds(1,4) = 1.06056257_wp; gauss_Ds(2,4) = 1.38282560_wp; gauss_Ds(3,4) = 2.40148179_wp; gauss_Ds(4,4) = 7.15513024_wp;

    real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
    gauss_wts(1,1) = 0.5_wp         ; gauss_wts(2,1) = 0._wp          ; gauss_wts(3,1) = 0._wp          ; gauss_wts(4,1) = 0._wp          ;
    gauss_wts(1,2) = 0.3180413817_wp; gauss_wts(2,2) = 0.1819586183_wp; gauss_wts(3,2) = 0._wp          ; gauss_wts(4,2) = 0._wp          ;
    gauss_wts(1,3) = 0.2009319137_wp; gauss_wts(2,3) = 0.2292411064_wp; gauss_wts(3,3) = 0.0698269799_wp; gauss_wts(4,3) = 0._wp          ;
    gauss_wts(1,4) = 0.1355069134_wp; gauss_wts(2,4) = 0.2034645680_wp; gauss_wts(3,4) = 0.1298475476_wp; gauss_wts(4,4) = 0.0311809710_wp;

    // Error checking
    //   if inc_flux is present it has the right dimensions, is positive definite
    int ncol  = optical_props.get_ncol();
    int nlay  = optical_props.get_nlay();
    int ngpt  = optical_props.get_ngpt();
    int nband = optical_props.get_nband();

    // Error checking -- consistency of sizes and validity of values
    if (! fluxes.are_desired()) { stoprun("rte_lw: no space allocated for fluxes"); }

    // Source functions
    if (sources.get_ncol() != ncol || sources.get_nlay() != lay || sources.get_ngpt() != ngpt) {
      stoprun("rte_lw: sources and optical properties inconsistently sized");
    }

    // Surface emissivity
    if (size(sfc_emis,1) != nband || size(sfc_emis,2) != ncol) { stoprun("rte_lw: sfc_emis inconsistently sized"); }
    if (anyLT(sfc_emis,0._wp) || anyGT(sfc_emis,1._wp)) { stoprun("rte_lw: sfc_emis has values < 0 or > 1"); }

    // Incident flux, if present
    if (allocated(inc_flux)) {
      if (size(inc_flux,1) != ncol | size(inc_flux,2) != ngpt) { stoprun("rte_lw: inc_flux inconsistently sized"); }
      if (anyLT(inc_flux,0._wp)) { stoprun("rte_lw: inc_flux has values < 0"); }
    }

    // Number of quadrature points for no-scattering calculation
    int n_quad_angs = 1;
    if ( n_gauss_angles != -1 ) {
      if (n_gauss_angles > max_gauss_pts) {
        stoprun("rte_lw: asking for too many quadrature points for no-scattering calculation");
      }
      if (n_gauss_angles < 1) {
        stoprun("rte_lw: have to ask for at least one quadrature point for no-scattering calculation");
      }
      n_quad_angs = n_gauss_angles;
    }
    
    // Ensure values of tau, ssa, and g are reasonable
    optical_props.validate();

    // Lower boundary condition -- expand surface emissivity by band to gpoints
    gpt_flux_up  = real3d("gpt_flux_up" ,ncol,nlay+1,ngpt);
    gpt_flux_dn  = real3d("gpt_flux_dn" ,ncol,nlay+1,ngpt);
    sfc_emis_gpt = real2d("sfc_emis_gpt",ncol       ,ngpt);
    expand_and_transpose(optical_props, sfc_emis, sfc_emis_gpt);
    
    //   Upper boundary condition
    if (allocated(inc_flux)) {
      apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, gpt_flux_dn);
    } else {
      // Default is zero incident diffuse flux
      apply_BC(ncol, nlay, ngpt, top_at_1          , gpt_flux_dn);
    }

    // Compute the radiative transfer...
    // No scattering two-stream calculation
    optical_props.validate();
    lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, logical(top_at_1, wl), &
                     n_quad_angs, gauss_Ds(1:n_quad_angs,n_quad_angs), gauss_wts(1:n_quad_angs,n_quad_angs), &
                     optical_props.tau,                                                  &
                     sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, &
                     sfc_emis_gpt, sources.sfc_source,  &
                     gpt_flux_up, gpt_flux_dn)
    // ...and reduce spectral fluxes to desired output quantities
    fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1)
  end function rte_lw


  // Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
  void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d &arr_out) {
    int2d limits;
    int ncol  = size(arr_in,2);
    int nband = ops.get_nband();
    int ngpt  = ops.get_ngpt();
    limits = ops.get_band_lims_gpoint();

    // for (int iband=1; iband <= nband; iband++) {
    //   for (int icol=1; icol <= ncol; icol++) {
    parallel_for( Bounds<2>(nband,ncol) , YAKL_LAMBDA (int iband, int icol) {
      for (int igpt=limits(1,iband); igpt <= limits(2,iband); igpt++) {
        arr_out(icol,igpt) = arr_in(iband,icol);
      }
    });
  }
}


