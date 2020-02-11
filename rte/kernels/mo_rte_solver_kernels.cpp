
#include "mo_rte_solver_kernels.h"



// -------------------------------------------------------------------------------------------------
//
//   Lower-level longwave kernels
//
// ---------------------------------------------------------------
//
// Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
// See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
// This routine implements point-wise stencil, and has to be called in a loop
//
// ---------------------------------------------------------------
YAKL_INLINE void lw_source_noscat_stencil(int ncol, int nlay, int ngpt, int icol, int ilay, int igpt,
                                          real3d const &lay_source, real3d const &lev_source_up, real3d const &lev_source_dn,
                                          real3d const &tau, real3d const &trans, real3d &source_dn, real3d &source_up) {
  real tau_thresh = sqrt( std::numeric_limits<real>::epsilon() );

  // ---------------------------------------------------------------
  //
  // Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
  //   is of order epsilon (smallest difference from 1. in working precision)
  //   Thanks to Peter Blossey
  //
  real fact = merge((1._wp - trans(icol,ilay,igpt))/tau(icol,ilay,igpt) - trans(icol,ilay,igpt), 
                        tau(icol,ilay,igpt) * ( 0.5_wp - 1._wp/3._wp*tau(icol,ilay,igpt) ), 
                        tau(icol,ilay,igpt) > tau_thresh);
  //
  // Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
  //
  source_dn(icol,ilay,igpt) = (1._wp - trans(icol,ilay,igpt)) * lev_source_dn(icol,ilay,igpt) + 
          2._wp * fact * (lay_source(icol,ilay,igpt) - lev_source_dn(icol,ilay,igpt));
  source_up(icol,ilay,igpt) = (1._wp - trans(icol,ilay,igpt)) * lev_source_up(icol,ilay,igpt) + 
          2._wp * fact * (lay_source(icol,ilay,igpt) - lev_source_up(icol,ilay,igpt));
}



// ---------------------
extern "C" void apply_BC_0(int ncol, int nlay, int ngpt, bool top_at_1, real *flux_dn_p) {
  umgReal3d flux_dn("flux_dn",flux_dn_p,ncol,nlay+1,ngpt);

  //   Upper boundary condition
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol,      1, igpt)  = 0;
    });
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol, nlay+1, igpt)  = 0;
    });
  }
}



extern "C" void apply_BC_factor(int ncol, int nlay, int ngpt, bool top_at_1, real *inc_flux_p, real *factor_p, real *flux_dn_p) {
  umgReal2d inc_flux("inc_flux",inc_flux_p,ncol       ,ngpt);
  umgReal2d factor  ("factor"  ,factor_p  ,ncol            );
  umgReal2d flux_dn ("flux_dn" ,flux_dn_p ,ncol,nlay+1,ngpt);

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt) * factor(icol);
    });
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt) * factor(icol);
    });
  }
}



// ---------------------------------------------------------------
//
// Upper boundary condition
//
// ---------------------------------------------------------------
extern "C" void apply_BC_gpt(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &inc_flux, real3d &flux_dn) {
  //   Upper boundary condition
  if (top_at_1) {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt);
    });
  } else {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt);
    });
  }
}



// ---------------------------------------------------------------
//
// Transport of diffuse radiation through a vertically layered atmosphere.
//   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
//   This routine is shared by longwave and shortwave
//
// -------------------------------------------------------------------------------------------------
extern "C" void adding(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &albedo_sfc,             
                       real3d const &rdif, real3d const &tdif, real3d const &src_dn, real3d const &src_up, real2d const &src_sfc, 
                       real3d &flux_up, real3d &flux_dn) {
  real3d albedo("albedo",ncol,nlay+1,ngpt);
  real3d src   ("src   ",ncol,nlay+1,ngpt);
  real3d denom ("denom ",ncol,nlay  ,ngpt);
  // ------------------
  // ---------------------------------
  //
  // Indexing into arrays for upward and downward propagation depends on the vertical
  //   orientation of the arrays (whether the domain top is at the first or last index)
  // We write the loops out explicitly so compilers will have no trouble optimizing them.
  //
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      int ilev = nlay + 1;
      // Albedo of lowest level is the surface albedo...
      albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
      // ... and source of diffuse radiation is surface emission
      src(icol,ilev,igpt) = src_sfc(icol,igpt);

      //
      // From bottom to top of atmosphere --
      //   compute albedo and source of upward radiation
      //
      for (ilev=nlay; ilev>=1; ilev--) {
        denom(icol,ilev,igpt) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(icol,ilev+1,igpt));    // Eq 10
        albedo(icol,ilev,igpt) = rdif(icol,ilev,igpt) + 
                                 tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev+1,igpt) * denom(icol,ilev,igpt); // Equation 9
        //
        // Equation 11 -- source is emitted upward radiation at top of layer plus
        //   radiation emitted at bottom of layer,
        //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
        //
        src(icol,ilev,igpt) =  src_up(icol, ilev, igpt) + 
                               tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *       
                               (src(icol,ilev+1,igpt) + albedo(icol,ilev+1,igpt)*src_dn(icol,ilev,igpt));
      }

      // Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = 1;
      flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                src(icol,ilev,igpt);                                  // emission from below

      //
      // From the top of the atmosphere downward -- compute fluxes
      //
      for (ilev = 1; ilev <= nlay+1; ilev++) {
        flux_dn(icol,ilev,igpt) = (tdif(icol,ilev-1,igpt)*flux_dn(icol,ilev-1,igpt) +   // Equation 13
                                  rdif(icol,ilev-1,igpt)*src(icol,ilev,igpt) +       
                                  src_dn(icol,ilev-1,igpt)) * denom(icol,ilev-1,igpt);
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                  src(icol,ilev,igpt);
      }
    });

  } else {

    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      int ilev = 1;
      // Albedo of lowest level is the surface albedo...
      albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
      // ... and source of diffuse radiation is surface emission
      src(icol,ilev,igpt) = src_sfc(icol,igpt);

      //
      // From bottom to top of atmosphere --
      //   compute albedo and source of upward radiation
      //
      for (ilev = 1; ilev <= nlay; ilev++) {
        denom (icol,ilev  ,igpt) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(icol,ilev,igpt));                // Eq 10
        albedo(icol,ilev+1,igpt) = rdif(icol,ilev,igpt) + 
                                   tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev,igpt) * denom(icol,ilev,igpt); // Equation 9
        //
        // Equation 11 -- source is emitted upward radiation at top of layer plus
        //   radiation emitted at bottom of layer,
        //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
        //
        src(icol,ilev+1,igpt) =  src_up(icol, ilev, igpt) +  
                                 tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *       
                                 (src(icol,ilev,igpt) + albedo(icol,ilev,igpt)*src_dn(icol,ilev,igpt));
      }

      // Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = nlay+1;
      flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                src(icol,ilev,igpt);                          // scattering by the direct beam below

      //
      // From the top of the atmosphere downward -- compute fluxes
      //
      for (ilev=nlay; ilev >= 1; ilev--) {
        flux_dn(icol,ilev,igpt) = (tdif(icol,ilev,igpt)*flux_dn(icol,ilev+1,igpt) +   // Equation 13
                                  rdif(icol,ilev,igpt)*src(icol,ilev,igpt) + 
                                  src_dn(icol, ilev, igpt)) * denom(icol,ilev,igpt);
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                  src(icol,ilev,igpt);

      }
    });
  }
}



// ---------------------------------------------------------------
//
// Direct beam source for diffuse radiation in layers and at surface;
//   report direct beam as a byproduct
//
YAKL_INLINE void sw_source_2str(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &Rdir, real3d const &Tdir,
                                real3d const &Tnoscat, real2d const &sfc_albedo, real3d &source_up, real3d &source_dn,
                                real2d &source_sfc, real3d &flux_dn_dir) {

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      for (int ilev=1; ilev<=nlay; ilev++) {
        source_up(icol,ilev,igpt)     =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        source_dn(icol,ilev,igpt)     =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        flux_dn_dir(icol,ilev+1,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        if (ilev == nlay) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,nlay+1,igpt)*sfc_albedo(icol,igpt);
        }
      }
    });
  } else {
    // layer index = level index
    // previous level is up (+1)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      for (int ilev=nlay; ilev>=1; ilev--) {
        source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        if (ilev ==    1) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,    1,igpt)*sfc_albedo(icol,igpt);
        }
      }
    });
  }
}



// -------------------------------------------------------------------------------------------------
//
//   Lower-level shortwave kernels
//
// -------------------------------------------------------------------------------------------------
//
// Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
//    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
//
// Equations are developed in Meador and Weaver, 1980,
//    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
//
YAKL_INLINE void sw_two_stream(int ncol, int nlay, int ngpt, real1d const &mu0, real3d const &tau,
                               real3d const &w0, real3d const &g, real3d &Rdif, real3d &Tdif,
                               real3d &Rdir, real3d &Tdir, real3d &Tnoscat) {
  real1d mu0_inv("mu0_inv",ncol);

  real eps = std::numeric_limits<real>::epsilon();

  parallel_for( Bounds<1>({1,ncol}) , YAKL_LAMBDA (int const indices[]) {
    int icol = indices[0];
    mu0_inv(icol) = 1._wp/mu0(icol);
  });

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );

    // Zdunkowski Practical Improved Flux Method "PIFM"
    //  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    //
    real gamma1= (8._wp - w0(icol,ilay,igpt) * (5._wp + 3._wp * g(icol,ilay,igpt))) * .25_wp;
    real gamma2=  3._wp *(w0(icol,ilay,igpt) * (1._wp -         g(icol,ilay,igpt))) * .25_wp;
    real gamma3= (2._wp - 3._wp * mu0(icol)  *                  g(icol,ilay,igpt) ) * .25_wp;
    real gamma4=  1._wp - gamma3;

    real alpha1 = gamma1 * gamma4 + gamma2 * gamma3;           // Eq. 16
    real alpha2 = gamma1 * gamma3 + gamma2 * gamma4;           // Eq. 17
    // Written to encourage vectorization of exponential, square root
    // Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
    //   k = 0 for isotropic, conservative scattering; this lower limit on k
    //   gives relative error with respect to conservative solution
    //   of < 0.1% in Rdif down to tau = 10^-9
    real k = sqrt(max((gamma1 - gamma2) * 
                      (gamma1 + gamma2),  
                      1.e-12_wp));
    real exp_minusktau = exp(-tau(icol,ilay,igpt)*k);
    //
    // Diffuse reflection and transmission
    //
    real exp_minus2ktau = exp_minusktau * exp_minusktau;

    // Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    real RT_term = 1._wp / (k      * (1._wp + exp_minus2ktau)  + 
                            gamma1 * (1._wp - exp_minus2ktau) );

    // Equation 25
    Rdif(icol,ilay,igpt) = RT_term * gamma2 * (1._wp - exp_minus2ktau);

    // Equation 26
    Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau;

    //
    // Transmittance of direct, unscattered beam. Also used below
    //
    Tnoscat(icol,ilay,igpt) = exp(-tau(icol,ilay,igpt)*mu0_inv(icol));

    //
    // Direct reflect and transmission
    //
    real k_mu     = k * mu0(icol);
    real k_gamma3 = k * gamma3;
    real k_gamma4 = k * gamma4;

    //
    // Equation 14, multiplying top and bottom by exp(-k*tau)
    //   and rearranging to avoid div by 0.
    //
    RT_term =  w0(icol,ilay,igpt) * RT_term/merge(1._wp - k_mu*k_mu, 
                                                  eps,    
                                                  abs(1._wp - k_mu*k_mu) >= eps);

    Rdir(icol,ilay,igpt) = RT_term  *                                    
       ((1._wp - k_mu) * (alpha2 + k_gamma3)                  - 
        (1._wp + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau - 
        2.0_wp * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau  * Tnoscat(icol,ilay,igpt));

    //
    // Equation 15, multiplying top and bottom by exp(-k*tau),
    //   multiplying through by exp(-tau/mu0) to
    //   prefer underflow to overflow
    // Omitting direct transmittance
    //
    Tdir(icol,ilay,igpt) = 
             -RT_term * ((1._wp + k_mu) * (alpha1 + k_gamma4) * Tnoscat(icol,ilay,igpt) - 
                         (1._wp - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat(icol,ilay,igpt) - 
                          2.0_wp * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau );

  });
}



// ---------------------------------------------------------------
//
// Longwave no-scattering transport
//
// ---------------------------------------------------------------
YAKL_INLINE void lw_transport_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &trans,
                                     real2d const &sfc_albedo, real3d const &source_dn, real3d const &source_up, real2d const &source_sfc, 
                                     real3d &radn_up, real3d &radn_dn) {
  if (top_at_1) {
    //
    // Top of domain is index 1
    //
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      // Downward propagation
      for (int ilev=2; ilev<=nlay+1; ilev++) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev-1,igpt)*radn_dn(icol,ilev-1,igpt) + source_dn(icol,ilev-1,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,nlay+1,igpt) = radn_dn(icol,nlay+1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=nlay; ilev>=1; ilev--) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_up(icol,ilev+1,igpt) + source_up(icol,ilev,igpt);
      }
    });
  } else {
    //
    // Top of domain is index nlay+1
    //
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      // Downward propagation
      for (int ilev=nlay; ilev>=1; ilev--) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,     1,igpt) = radn_dn(icol,     1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=2; ilev<=nlay+1; ilev++) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev-1,igpt) * radn_up(icol,ilev-1,igpt) +  source_up(icol,ilev-1,igpt);
      }
    });
  }
}



extern "C" void sw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, real *tau_p, real *ssa_p, real *g_p,
                                  real *mu0_p, real *sfc_alb_dir_p, real *sfc_alb_dif_p, real *flux_up_p,
                                  real *flux_dn_p, real *flux_dir_p) {
  umgReal3d tau         ("tau        ",tau_p        ,ncol,nlay,  ngpt);
  umgReal3d ssa         ("ssa        ",ssa_p        ,ncol,nlay,  ngpt);
  umgReal3d g           ("g          ",g_p          ,ncol,nlay,  ngpt);
  umgReal1d mu0         ("mu0        ",mu0_p        ,ncol            );
  umgReal2d sfc_alb_dir ("sfc_alb_dir",sfc_alb_dir_p,ncol,       ngpt);
  umgReal2d sfc_alb_dif ("sfc_alb_dif",sfc_alb_dif_p,ncol,       ngpt);
  umgReal3d flux_up     ("flux_up    ",flux_up_p    ,ncol,nlay+1,ngpt);
  umgReal3d flux_dn     ("flux_dn    ",flux_dn_p    ,ncol,nlay+1,ngpt);
  umgReal3d flux_dir    ("flux_dir   ",flux_dir_p   ,ncol,nlay+1,ngpt);
  // -------------------------------------------
  real3d Rdif      ("Rdif      ",ncol,nlay,ngpt);         
  real3d Tdif      ("Tdif      ",ncol,nlay,ngpt);         
  real3d Rdir      ("Rdir      ",ncol,nlay,ngpt);         
  real3d Tdir      ("Tdir      ",ncol,nlay,ngpt);         
  real3d Tnoscat   ("Tnoscat   ",ncol,nlay,ngpt);         
  real3d source_up ("source_up ",ncol,nlay,ngpt);         
  real3d source_dn ("source_dn ",ncol,nlay,ngpt);         
  real2d source_srf("source_srf",ncol     ,ngpt);         
  // ------------------------------------
  //
  // Cell properties: transmittance and reflectance for direct and diffuse radiation
  //
  sw_two_stream(ncol, nlay, ngpt, mu0, 
                tau , ssa , g   ,      
                Rdif, Tdif, Rdir, Tdir, Tnoscat);

  sw_source_2str(ncol, nlay, ngpt, top_at_1,       
                 Rdir, Tdir, Tnoscat, sfc_alb_dir, 
                 source_up, source_dn, source_srf, flux_dir);

  adding(ncol, nlay, ngpt, top_at_1,   
         sfc_alb_dif, Rdif, Tdif,      
         source_dn, source_up, source_srf, flux_up, flux_dn);

  //
  // adding computes only diffuse flux; flux_dn is total
  //
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay+1
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay+1},{1,ncol}) , YAKL_LAMBDA (int const indices[]) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );
    flux_dn(icol,ilay,igpt) = flux_dn(icol,ilay,igpt) + flux_dir(icol,ilay,igpt);
  });
}



// -------------------------------------------------------------------------------------------------
//
// Top-level longwave kernels
//
// -------------------------------------------------------------------------------------------------
//
// LW fluxes, no scattering, mu (cosine of integration angle) specified by column
//   Does radiation calculation at user-supplied angles; converts radiances to flux
//   using user-supplied weights
//
// ---------------------------------------------------------------
extern "C" void lw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &D, real weight,
                                 real3d const &tau, real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                                 real2d const &sfc_emis, real2d const &sfc_src, real3d &radn_up, real3d &radn_dn) {
  real3d tau_loc   ("tau_loc   ",ncol,nlay,ngpt);             
  real3d trans     ("trans     ",ncol,nlay,ngpt);             
  real3d source_dn ("source_dn ",ncol,nlay,ngpt);             
  real3d source_up ("source_up ",ncol,nlay,ngpt);             
  real2d source_sfc("source_sfc",ncol,     ngpt);             
  real2d sfc_albedo("sfc_albedo",ncol,     ngpt);             

  real pi = acos(-1._wp);

  // ------------------------------------
  // Which way is up?
  // Level Planck sources for upward and downward radiation
  // When top_at_1, lev_source_up => lev_source_dec
  //                lev_source_dn => lev_source_inc, and vice-versa
  int top_level;
  umgReal3d lev_source_up;
  umgReal3d lev_source_dn;
  if (top_at_1) {
    top_level = 1;
    lev_source_up = umgReal3d("lev_source_up",lev_source_dec.data(),ncol,nlay,ngpt);
    lev_source_dn = umgReal3d("lev_source_dn",lev_source_inc.data(),ncol,nlay,ngpt);
  } else {
    top_level = nlay+1;
    lev_source_up = umgReal3d("lev_source_up",lev_source_inc.data(),ncol,nlay,ngpt);
    lev_source_dn = umgReal3d("lev_source_dn",lev_source_dec.data(),ncol,nlay,ngpt);
  }

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
    int igpt, icol;
    storeIndices( indices , igpt,icol );

    // Transport is for intensity
    //   convert flux at top of domain to intensity assuming azimuthal isotropy
    //
    radn_dn(icol,top_level,igpt) = radn_dn(icol,top_level,igpt)/(2._wp * pi * weight);
    
    // Surface albedo, surface source function
    //
    sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt);
    source_sfc(icol,igpt) = sfc_emis(icol,igpt) * sfc_src(icol,igpt);
  });

  // NOTE: This kernel produces small differences between GPU and CPU
  // implementations on Ascent with PGI, we assume due to floating point
  // differences in the exp() function. These differences are small in the
  // RFMIP test case (10^-6).
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay},{1,ncol}) , YAKL_LAMBDA (int const indices[]) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );
    //
    // Optical path and transmission, used in source function and transport calculations
    //
    tau_loc(icol,ilay,igpt) = tau(icol,ilay,igpt)*D(icol,igpt);
    trans  (icol,ilay,igpt) = exp(-tau_loc(icol,ilay,igpt));

    lw_source_noscat_stencil(ncol, nlay, ngpt, icol, ilay, igpt,        
                             lay_source, lev_source_up, lev_source_dn,  
                             tau_loc, trans,                            
                             source_dn, source_up);
  });

  //
  // Transport
  //
  lw_transport_noscat(ncol, nlay, ngpt, top_at_1,  
                      tau_loc, trans, sfc_albedo, source_dn, source_up, source_sfc, 
                      radn_up, radn_dn);

  //
  // Convert intensity to flux assuming azimuthal isotropy and quadrature weight
  //
  //$acc parallel loop collapse(3)
  // do igpt = 1, ngpt
  //   do ilev = 1, nlay+1
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay+1},{1,ncol}) , YAKL_LAMBDA (int const indices[]) {
    int igpt, ilev, icol;
    storeIndices( indices , igpt,ilev,icol );

    radn_dn(icol,ilev,igpt) = 2._wp * pi * weight * radn_dn(icol,ilev,igpt);
    radn_up(icol,ilev,igpt) = 2._wp * pi * weight * radn_up(icol,ilev,igpt);
  });
}



// ---------------------------------------------------------------
//
// LW transport, no scattering, multi-angle quadrature
//   Users provide a set of weights and quadrature angles
//   Routine sums over single-angle solutions for each sets of angles/weights
//
// ---------------------------------------------------------------
extern "C" void lw_solver_noscat_GaussQuad(int ncol, int nlay, int ngpt, bool top_at_1, int nmus, real *Ds_p, real *weights_p, 
                                           real *tau_p, real *lay_source_p, real *lev_source_inc_p, real *lev_source_dec_p,
                                           real *sfc_emis_p, real *sfc_src_p, real *flux_up_p, real *flux_dn_p) {

  umgReal1d Ds             ("Ds            ",Ds_p            ,nmus            );
  umgReal1d weights        ("weights       ",weights_p       ,nmus            );
  umgReal3d tau            ("tau           ",tau_p           ,ncol,nlay,  ngpt);
  umgReal3d lay_source     ("lay_source    ",lay_source_p    ,ncol,nlay,  ngpt);
  umgReal3d lev_source_inc ("lev_source_inc",lev_source_inc_p,ncol,nlay,  ngpt);
  umgReal3d lev_source_dec ("lev_source_dec",lev_source_dec_p,ncol,nlay,  ngpt);
  umgReal2d sfc_emis       ("sfc_emis      ",sfc_emis_p      ,ncol,       ngpt);
  umgReal2d sfc_src        ("sfc_src       ",sfc_src_p       ,ncol,       ngpt);
  umgReal3d flux_dn        ("flux_dn       ",flux_dn_p       ,ncol,nlay+1,ngpt);
  umgReal3d flux_up        ("flux_up       ",flux_up_p       ,ncol,nlay+1,ngpt);

  // Local variables
  real3d radn_dn ("radn_dn ",ncol,nlay+1,ngpt);  
  real3d radn_up ("radn_up ",ncol,nlay+1,ngpt);  
  real2d Ds_ncol ("Ds_ncol ",ncol,       ngpt);  
  real2d flux_top("flux_top",ncol,       ngpt);  

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
    int igpt, icol;
    storeIndices( indices , igpt,icol );

    Ds_ncol(icol, igpt) = Ds(1);
  });

  lw_solver_noscat(ncol, nlay, ngpt, 
                   top_at_1, Ds_ncol, weights(1), tau, 
                   lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, 
                   flux_up, flux_dn);
  //
  // For more than one angle use local arrays
  //
  int top_level = merge(1, nlay+1, top_at_1);

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
    int igpt, icol;
    storeIndices( indices , igpt,icol );

    flux_top(icol,igpt) = flux_dn(icol,top_level,igpt);
  });

  apply_BC_gpt(ncol, nlay, ngpt, top_at_1, flux_top, radn_dn);

  for (int imu=2; imu<=nmus; imu++) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, icol;
      storeIndices( indices , igpt,icol );

      Ds_ncol(icol, igpt) = Ds(imu);
    });

    lw_solver_noscat(ncol, nlay, ngpt, 
                     top_at_1, Ds_ncol, weights(imu), tau, 
                     lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, 
                     radn_up, radn_dn);

    // do igpt = 1, ngpt
    //   do ilev = 1, nlay+1
    //     do icol = 1, ncol
    parallel_for( Bounds<3>({1,ngpt},{1,nlay+1},{1,ncol}) , YAKL_LAMBDA ( int const indices[] ) {
      int igpt, ilev, icol;
      storeIndices( indices , igpt,ilev,icol );

      flux_up(icol,ilev,ngpt) = flux_up(icol,ilev,ngpt) + radn_up(icol,ilev,ngpt);
      flux_dn(icol,ilev,ngpt) = flux_dn(icol,ilev,ngpt) + radn_dn(icol,ilev,ngpt);
    });

  } // imu
}




