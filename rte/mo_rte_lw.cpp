
#include "mo_rte_lw.h"

// Interface using only optical properties and source functions as inputs; fluxes as outputs.
void rte_lw(OpticalProps1scl const &optical_props, bool top_at_1, SourceFuncLW const &sources, real2d const &sfc_emis,
            FluxesBroadband &fluxes, real2d const &inc_flux, int n_gauss_angles) {
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
  if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
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
  real1d tmp_Ds ("tmp_Ds" ,n_quad_angs);
  real1d tmp_wts("tmp_wts",n_quad_angs);
  // for (int i=1 ; i <= n_quad_angs ; i++) {
  parallel_for( Bounds<1>(n_quad_angs) , YAKL_LAMBDA (int i) {
    tmp_Ds (i) = gauss_Ds (i,n_quad_angs);
    tmp_wts(i) = gauss_wts(i,n_quad_angs);
  });
  lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, n_quad_angs, tmp_Ds, tmp_wts, optical_props.tau,                                                  
                             sources.lay_source, sources.lev_source_inc, sources.lev_source_dec, 
                             sfc_emis_gpt, sources.sfc_source, gpt_flux_up, gpt_flux_dn);
  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);
}


