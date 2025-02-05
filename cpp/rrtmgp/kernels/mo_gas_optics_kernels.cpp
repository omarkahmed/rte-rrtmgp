
#include "mo_gas_optics_kernels.h"
#include <limits>

using yakl::SB;
using yakl::COLON;


void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, int2d const &flavor,
                   real1d const &press_ref_log, real1d const &temp_ref, real press_ref_log_delta, real temp_ref_min,
                   real temp_ref_delta, real press_ref_trop_log, real3d const &vmr_ref, real2d const &play,
                   real2d const &tlay, real3d const &col_gas, int2d &jtemp, real6d &fmajor, real5d &fminor,
                   real4d &col_mix, bool2d &tropo, int4d &jeta, int2d &jpress) {

  real2d ftemp ("ftemp" ,ncol,nlay);
  real2d fpress("fpress",ncol,nlay);

  real tiny = std::numeric_limits<real>::min();

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  parallel_for( Bounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
    // index and factor for temperature interpolation
    jtemp(icol,ilay) = (int) ((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
    jtemp(icol,ilay) = min(ntemp - 1, max(1, jtemp(icol,ilay))); // limit the index range
    ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta;

    // index and factor for pressure interpolation
    real locpress = 1. + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta;
    jpress(icol,ilay) = min(npres-1, max(1, (int)(locpress)));
    fpress(icol,ilay) = locpress - (real)(jpress(icol,ilay));

    // determine if in lower or upper part of atmosphere
    tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log;
  });

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int iflav=1; iflav<=nflav; iflav++) {   // loop over implemented combinations of major species
  //       for (int itemp=1; itemp<=2; itemp++) {
  parallel_for( Bounds<4>(nlay,ncol,nflav,2) , YAKL_LAMBDA (int ilay, int icol, int iflav , int itemp) {
    yakl::FSArray<int,1,SB<2>> igases;

    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));
    igases(1) = flavor(1,iflav);
    igases(2) = flavor(2,iflav);

    // compute interpolation fractions needed for lower, then upper reference temperature level
    // compute binary species parameter (eta) for flavor and temperature and
    //  associated interpolation index and factors
    real ratio_eta_half = vmr_ref(itropo,igases(1),(jtemp(icol,ilay)+itemp-1)) / 
                          vmr_ref(itropo,igases(2),(jtemp(icol,ilay)+itemp-1));
    col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2));
    real eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(itemp,iflav,icol,ilay), 0.5, 
                     col_mix(itemp,iflav,icol,ilay) > 2. * tiny);
    real loceta = eta * (neta-1.0);
    jeta(itemp,iflav,icol,ilay) = min((int)(loceta)+1, neta-1);
    real feta = fmod(loceta, 1.0);
    // compute interpolation fractions needed for minor species
    real ftemp_term = ((2.0 - itemp) + (2.0 * itemp - 3.0 ) * ftemp(icol,ilay));
    fminor(1,itemp,iflav,icol,ilay) = (1. - feta) * ftemp_term;
    fminor(2,itemp,iflav,icol,ilay) =          feta  * ftemp_term;
    // compute interpolation fractions needed for major species
    fmajor(1,1,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,1,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(2,itemp,iflav,icol,ilay);
    fmajor(1,2,itemp,iflav,icol,ilay) =          fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,2,itemp,iflav,icol,ilay) =          fpress(icol,ilay)  * fminor(2,itemp,iflav,icol,ilay);
  });
}



void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real3d const &tau_abs, real3d const &tau_rayleigh,
                              real3d &tau, real3d &ssa, real3d &g) {

  real tiny = std::numeric_limits<real>::min();

  int constexpr TILE_SIZE=2;
  int colTiles = ncol / TILE_SIZE + 1;
  int gptTiles = ngpt / TILE_SIZE + 1;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int tcol=1; tcol<=colTiles; tcol++) {
  //     for (int tgpt=1; tgpt<=gptTiles; tgpt++) {
  //       for (int itcol=1; itcol<=TILE_SIZE; itcol++) {
  //         for (int itgpt=1; itgpt<=TILE_SIZE; itgpt++) {
  parallel_for( Bounds<5>(nlay,colTiles,gptTiles,TILE_SIZE,TILE_SIZE) , YAKL_LAMBDA (int ilay, int tcol, int tgpt, int itcol, int itgpt) {
    int icol = (tcol-1)*TILE_SIZE + itcol;
    int igpt = (tgpt-1)*TILE_SIZE + itgpt;

    if ( icol <= ncol && igpt <= ngpt ) {
      real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
      tau(icol,ilay,igpt) = t;
      g  (icol,ilay,igpt) = 0.;
      if(t > 2. * tiny) {
        ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
      } else {
        ssa(icol,ilay,igpt) = 0.;
      }
    }
  });
}



void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           real2d const &tlay, real2d const &tlev, real1d const &tsfc, int sfc_lay, real6d const &fmajor,
                           int4d const &jeta, bool2d const &tropo, int2d const &jtemp, int2d const &jpress,
                           int1d const &gpoint_bands, int2d const &band_lims_gpt, real4d const &pfracin, real temp_ref_min,
                           real totplnk_delta, real2d const &totplnk, int2d const &gpoint_flavor, real2d &sfc_src,
                           real3d &lay_src, real3d &lev_src_inc, real3d &lev_src_dec) {

  real3d pfrac          ("pfrac"          ,ngpt,nlay,ncol);
  real3d planck_function("planck_function",nbnd,nlay+1,ncol);
  real1d one            ("one"            ,2);
  memset(one,1.);

  // Calculation of fraction of band's Planck irradiance associated with each g-point
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( Bounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
    int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
    // interpolation in temperature, pressure, and eta
    pfrac(igpt,ilay,icol) = 
      interpolate3D(one, fmajor.slice<3>(COLON,COLON,COLON,iflav,icol,ilay), pfracin, 
                    igpt, jeta.slice<1>(COLON,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp);
  });

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  // for (int icol=1; icol<=ncol; icol++) {
  parallel_for( ncol , YAKL_LAMBDA (int icol) {
    real1d planck_function_slice = planck_function.slice<1>(COLON,1,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // for (int igpt=1; igpt<=ngpt; igpt++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  parallel_for( Bounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
    sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 1, icol);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  parallel_for( Bounds<2>(ncol,nlay) , YAKL_LAMBDA (int icol, int ilay) {
    // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
    real1d planck_function_slice = planck_function.slice<1>(COLON,ilay,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // Explicitly unroll a time-consuming loop here to increase instruction-level parallelism on a GPU
  // Helps to achieve higher bandwidth
  //
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( Bounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
  });

  // compute level source irradiances for each g-point, one each for upward and downward paths
  // for (int icol=1; icol<=ncol; icol++) {
  parallel_for( ncol , YAKL_LAMBDA (int icol) {
    real1d planck_function_slice = planck_function.slice<1>(COLON,1,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tlev(icol,1), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=2; ilay<=nlay+1; ilay++) {
  parallel_for( Bounds<2>(ncol,{2,nlay+1}) , YAKL_LAMBDA (int icol, int ilay) {
    real1d planck_function_slice = planck_function.slice<1>(COLON,ilay,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });

  //
  // Map to g-points
  //
  // Same unrolling as mentioned before
  //
  // for (int icol=1; icol<=ncol; icol+=2) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( Bounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
    lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
    if (icol < ncol) {
      lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
      lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
    }
  });
}



// compute Rayleigh scattering optical depths
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &krayl, int idx_h2o,
                          real2d const &col_dry, real3d const &col_gas, real5d const &fminor, int4d const &jeta,
                          bool2d const &tropo, int2d const &jtemp, real3d &tau_rayleigh) {

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( Bounds<3>(nlay,ncol,ngpt) , YAKL_LAMBDA (int ilay, int icol, int igpt) {
    int itropo = merge(1,2,tropo(icol,ilay)); // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int iflav = gpoint_flavor(itropo, igpt);
    real k = interpolate2D(fminor.slice<2>(COLON,COLON,iflav,icol,ilay), 
                           krayl.slice<3>(COLON,COLON,COLON,itropo),      
                           igpt, jeta.slice<1>(COLON,iflav,icol,ilay), jtemp(icol,ilay), ngpt, neta, ntemp);
    tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay));
  });
}



// compute minor species optical depths
void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2d const &gpt_flv,
                              real3d const &kminor, int2d const &minor_limits_gpt, bool1d const &minor_scales_with_density,
                              bool1d const &scale_by_complement, int1d const &idx_minor, int1d const &idx_minor_scaling,
                              int1d const &kminor_start, real2d const &play, real2d const &tlay, real3d const &col_gas,
                              real5d const &fminor, int4d const &jeta, int2d const &layer_limits, int2d const &jtemp, real3d &tau) {
  real constexpr PaTohPa = 0.01;

  int extent = nminor;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt0=0; igpt0<=max_gpt_diff; igpt0++) {
  parallel_for( Bounds<3>(nlay,ncol,{0,max_gpt_diff}) , YAKL_DEVICE_LAMBDA (int ilay, int icol, int igpt0) {
    // This check skips individual columns with no pressures in range
    //
    if ( layer_limits(icol,1) <= 0 || ilay < layer_limits(icol,1) || ilay > layer_limits(icol,2) ) {
    } else {
      real myplay  = play (icol,ilay);
      real mytlay  = tlay (icol,ilay);
      real myjtemp = jtemp(icol,ilay);
      real mycol_gas_h2o = col_gas(icol,ilay,idx_h2o);
      real mycol_gas_0   = col_gas(icol,ilay,0);

      for (int imnr=1; imnr<=extent; imnr++) {

        real scaling = col_gas(icol,ilay,idx_minor(imnr));
        if (minor_scales_with_density(imnr)) {
          //
          // NOTE: P needed in hPa to properly handle density scaling.
          //
          scaling = scaling * (PaTohPa * myplay/mytlay);

          if (idx_minor_scaling(imnr) > 0) {  // there is a second gas that affects this gas's absorption
            real mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr));
            real vmr_fact = 1. / mycol_gas_0;
            real dry_fact = 1. / (1. + mycol_gas_h2o * vmr_fact);
            // scale by density of special gas
            if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
              scaling = scaling * (1. - mycol_gas_imnr * vmr_fact * dry_fact);
            } else {
              scaling = scaling *          mycol_gas_imnr * vmr_fact * dry_fact;
            }
          }
        } // minor_scalse_with_density(imnr)

        //
        // Interpolation of absorption coefficient and calculation of optical depth
        //
        // Which gpoint range does this minor gas affect?
        int gptS = minor_limits_gpt(1,imnr);
        int gptE = minor_limits_gpt(2,imnr);

        // Find the actual g-point to work on
        int igpt = igpt0 + gptS;

        // Proceed only if the g-point is within the correct range
        if (igpt <= gptE) {
          // What is the starting point in the stored array of minor absorption coefficients?
          int minor_start = kminor_start(imnr);

          real tau_minor = 0.;
          int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
          int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
          real kminor_loc = interpolate2D(fminor.slice<2>(COLON,COLON,iflav,icol,ilay), kminor, minor_loc,  
                                          jeta.slice<1>(COLON,iflav,icol,ilay), myjtemp, nminork, neta, ntemp);
          tau_minor = kminor_loc * scaling;

          yakl::atomicAdd( tau(igpt,ilay,icol) , tau_minor );
        }  // igpt <= gptE
      }
    }
  });
}



// compute minor species optical depths
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &kmajor,                         
                              real4d const &col_mix, real6d const &fmajor, int4d const &jeta, bool2d const &tropo, 
                              int2d const &jtemp, int2d const &jpress, real3d &tau) {
  // optical depth calculation for major species
  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     // optical depth calculation for major species
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  parallel_for( Bounds<3>(nlay,ncol,ngpt) , YAKL_LAMBDA (int ilay, int icol, int igpt) {
    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  // WS: moved inside innermost loop

    // binary species parameter (eta) and col_mix depend on band flavor
    int iflav = gpoint_flavor(itropo, igpt);
    real tau_major = 
      // interpolation in temperature, pressure, and eta
      interpolate3D(col_mix.slice<1>(COLON,iflav,icol,ilay), 
                    fmajor.slice<3>(COLON,COLON,COLON,iflav,icol,ilay), kmajor, 
                    igpt, jeta.slice<1>(COLON,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp);
    tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_major;
  });
}



// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, int2d const &gpoint_flavor,
                            int2d const &band_lims_gpt, real4d const &kmajor, real3d const &kminor_lower, real3d const &kminor_upper,
                            int2d const &minor_limits_gpt_lower, int2d const &minor_limits_gpt_upper, bool1d const &minor_scales_with_density_lower,    
                            bool1d const &minor_scales_with_density_upper, bool1d const &scale_by_complement_lower,          
                            bool1d const &scale_by_complement_upper, int1d const &idx_minor_lower, int1d const &idx_minor_upper,
                            int1d const &idx_minor_scaling_lower, int1d const &idx_minor_scaling_upper, int1d const &kminor_start_lower,                 
                            int1d const &kminor_start_upper, bool2d const &tropo, real4d const &col_mix, real6d const &fmajor,
                            real5d const &fminor, real2d const &play, real2d const &tlay, real3d const &col_gas, int4d const &jeta,
                            int2d const &jtemp, int2d const &jpress, real3d &tau, bool top_at_1) {

  int2d itropo_lower("itropo_lower",ncol,2);
  int2d itropo_upper("itropo_upper",ncol,2);

  int huge  = std::numeric_limits<int>::max();
  int small = std::numeric_limits<int>::min();

  if (top_at_1) {

    // for (int icol=1; icol<=ncol; icol++){
    parallel_for( ncol , YAKL_LAMBDA (int icol) {
      itropo_lower(icol,2) = nlay;
      // itropo_lower(icol,1) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=1; i<=nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,1) = minloc;
      }

      itropo_upper(icol,1) = 1;
      // itropo_upper(icol,2) = maxloc(play(icol,:), dim=1, mask=(.not. tropo(icol,:)))
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=1; i<=nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,2) = maxloc;
      }
    });

  } else {  // top_at_1

    // for (int icol=1; icol<=ncol; icol++){
    parallel_for( ncol , YAKL_LAMBDA ( int icol ) {
      itropo_lower(icol,1) = 1;
      // itropo_lower(icol,2) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
      {
        int minloc = huge;
        real mn = (real) huge;
        for (int i=1; i<=nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = min(mn,play(icol,i));
          }
        }
        itropo_lower(icol,2) = minloc;
      }

      itropo_upper(icol,2) = nlay;
      // itropo_upper(icol,1) = maxloc(play(icol,:), dim=1, mask=(.not.tropo(icol,:)))
      {
        int maxloc = small;
        real mx = (real) small;
        for (int i=1; i<=nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = max(mx,play(icol,i));
          }
        }
        itropo_upper(icol,1) = maxloc;
      }

    });

  }  // top_at_1

  // ---------------------
  // Major Species
  // ---------------------
  gas_optical_depths_major(ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, gpoint_flavor,
                           band_lims_gpt, kmajor, col_mix, fmajor, jeta, tropo, jtemp, jpress, tau);
  // ---------------------
  // Minor Species - lower
  // ---------------------
  int idx_tropo = 1;
  gas_optical_depths_minor(max_gpt_diff_lower, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorlower, nminorklower,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_lower, minor_limits_gpt_lower,
                           minor_scales_with_density_lower, scale_by_complement_lower, idx_minor_lower,
                           idx_minor_scaling_lower, kminor_start_lower, play,  tlay, col_gas, fminor,
                           jeta, itropo_lower, jtemp, tau);
  // ---------------------
  // Minor Species - upper
  // ---------------------
  idx_tropo = 2;
  gas_optical_depths_minor(max_gpt_diff_upper, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorupper, nminorkupper,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_upper, minor_limits_gpt_upper,
                           minor_scales_with_density_upper, scale_by_complement_upper, idx_minor_upper,
                           idx_minor_scaling_upper, kminor_start_upper, play, tlay, col_gas, fminor,
                           jeta, itropo_upper, jtemp, tau);
}



// Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
//   using Rayleigh scattering phase function
void combine_and_reorder_nstr(int ncol, int nlay, int ngpt, int nmom, real3d &tau_abs, real3d &tau_rayleigh, real3d &tau, real3d &ssa, real4d &p) {
  real tiny = std::numeric_limits<real>::min();

  // do icol = 1, ncol
  //   do ilay = 1, nlay
  //     do igpt = 1, ngpt
  parallel_for( Bounds<3>(ncol,nlay,ngpt) , YAKL_LAMBDA (int icol, int ilay, int igpt) {
    real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
    tau(icol,ilay,igpt) = t;
    if (t > 2. * tiny) {
      ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
    } else {
      ssa(icol,ilay,igpt) = 0;
    }
    for (int imom=1; imom<=nmom; imom++) {
      p(icol,ilay,igpt,imom) = 0;
    }
    if (nmom >= 2) {
      p(icol,ilay,igpt,2) = 0.1;
    }
  });
}



