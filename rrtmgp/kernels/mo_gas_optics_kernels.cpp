
#include "mo_gas_optics_kernels.h"
#include <limits>

using yakl::bnd;
using yakl::COLON;
using yakl::Bounds;


extern "C" void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp,
                              int *flavor_p, real *press_ref_log_p, real *temp_ref_p, real press_ref_log_delta,
                              real temp_ref_min, real temp_ref_delta, real press_ref_trop_log,
                              real *vmr_ref_p, real *play_p, real *tlay_p, real *col_gas_p, int *jtemp_p,
                              real *fmajor_p, real *fminor_p, real *col_mix_p, bool *tropo_p, int *jeta_p,
                              int *jpress_p) {

  umgInt2d  flavor       ("flavor"       ,flavor_p       ,2,nflav);
  umgReal1d press_ref_log("press_ref_log",press_ref_log_p,npres);
  umgReal1d temp_ref     ("temp_ref"     ,temp_ref_p     ,ntemp);
  umgReal3d vmr_ref      ("vmr_ref"      ,vmr_ref_p      , {1,2} , {0,ngas} , {1,ntemp} );
  umgReal2d play         ("play"         ,play_p         ,ncol,nlay);
  umgReal2d tlay         ("tlay"         ,tlay_p         ,ncol,nlay);
  umgReal3d col_gas      ("col_gas"      ,col_gas_p      , {1,ncol} , {1,nlay} , {0,ngas} );
  umgInt2d  jtemp        ("jtemp"        ,jtemp_p        ,ncol,nlay);
  umgInt2d  jpress       ("jpress"       ,jpress_p       ,ncol,nlay);
  umgBool2d tropo        ("tropo"        ,tropo_p        ,ncol,nlay);
  umgInt4d  jeta         ("jeta"         ,jeta_p         ,2,    nflav,ncol,nlay);
  umgReal4d col_mix      ("col_mix"      ,col_mix_p      ,2,    nflav,ncol,nlay);
  umgReal6d fmajor       ("fmajor"       ,fmajor_p       ,2,2,2,nflav,ncol,nlay);
  umgReal5d fminor       ("fminor"       ,fminor_p       ,2,2,  nflav,ncol,nlay);

  real2d ftemp ("ftemp" ,ncol,nlay);
  real2d fpress("fpress",ncol,nlay);

  real tiny = std::numeric_limits<real>::min();

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  yakl::parallel_for_cpu_serial( Bounds<2>({1,nlay},{1,ncol}) , YAKL_LAMBDA (int indices[]) {
    int ilay, icol;
    yakl::storeIndices( indices , ilay,icol );

    // index and factor for temperature interpolation
    jtemp(icol,ilay) = (int) ((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
    jtemp(icol,ilay) = min(ntemp - 1, max(1, jtemp(icol,ilay))); // limit the index range
    ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta;

    // index and factor for pressure interpolation
    real locpress = 1._wp + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta;
    jpress(icol,ilay) = min(npres-1, max(1, (int)(locpress)));
    fpress(icol,ilay) = locpress - (real)(jpress(icol,ilay));

    // determine if in lower or upper part of atmosphere
    tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log;
  });

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int iflav=1; iflav<=nflav; iflav++) {   // loop over implemented combinations of major species
  //       for (int itemp=1; itemp<=2; itemp++) {
  yakl::parallel_for_cpu_serial( Bounds<4>({1,nlay},{1,ncol},{1,nflav},{1,2}) , YAKL_LAMBDA ( int indices[]) {
    int ilay, icol, iflav, itemp;
    yakl::storeIndices( indices , ilay,icol,iflav,itemp );

    yakl::FSArray<int,bnd<1,2>> igases;

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
    real eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(itemp,iflav,icol,ilay), 0.5_wp, 
                     col_mix(itemp,iflav,icol,ilay) > 2._wp * tiny);
    real loceta = eta * (neta-1.0_wp);
    jeta(itemp,iflav,icol,ilay) = min((int)(loceta)+1, neta-1);
    real feta = fmod(loceta, 1.0_wp);
    // compute interpolation fractions needed for minor species
    real ftemp_term = ((2.0_wp-itemp) + (2.0_wp*itemp-3.0_wp) * ftemp(icol,ilay));
    fminor(1,itemp,iflav,icol,ilay) = (1._wp-feta) * ftemp_term;
    fminor(2,itemp,iflav,icol,ilay) =        feta  * ftemp_term;
    // compute interpolation fractions needed for major species
    fmajor(1,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(2,itemp,iflav,icol,ilay);
    fmajor(1,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay);
    fmajor(2,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(2,itemp,iflav,icol,ilay);
  });
}



extern "C" void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real *tau_abs_p, real *tau_rayleigh_p,
                                         real *tau_p,  real *ssa_p, real *g_p) {

  umgReal3d tau_abs     ("tau_abs"     ,tau_abs_p     ,ngpt,nlay,ncol);
  umgReal3d tau_rayleigh("tau_rayleigh",tau_rayleigh_p,ngpt,nlay,ncol);
  umgReal3d tau         ("tau"         ,tau_p         ,ncol,nlay,ngpt);
  umgReal3d ssa         ("ssa"         ,ssa_p         ,ncol,nlay,ngpt);
  umgReal3d g           ("g"           ,g_p           ,ncol,nlay,ngpt);

  real tiny = std::numeric_limits<real>::min();

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  yakl::parallel_for_cpu_serial( Bounds<3>({1,ncol},{1,nlay},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int icol,ilay,igpt;
    yakl::storeIndices( indices , icol,ilay,igpt );
     real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
     tau(icol,ilay,igpt) = t;
     g  (icol,ilay,igpt) = 0._wp;
     if(t > 2._wp * tiny) {
       ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
     } else {
       ssa(icol,ilay,igpt) = 0._wp;
     }
  });
}



extern "C" void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta,
                                      int npres, int ntemp, int nPlanckTemp, real *tlay_p, real *tlev_p,
                                      real *tsfc_p, int sfc_lay, real *fmajor_p, int *jeta_p, bool *tropo_p,
                                      int *jtemp_p, int *jpress_p, int *gpoint_bands_p, int *band_lims_gpt_p,           
                                      real *pfracin_p, real temp_ref_min, real totplnk_delta, real *totplnk_p,
                                      int *gpoint_flavor_p, real *sfc_src_p, real *lay_src_p, real *lev_src_inc_p,
                                      real *lev_src_dec_p) {

  umgReal2d tlay          ("tlay"          ,tlay_p          ,ncol,nlay);
  umgReal2d tlev          ("tlev"          ,tlev_p          ,ncol,nlay+1);
  umgReal1d tsfc          ("tsfc"          ,tsfc_p          ,ncol);
  umgReal6d fmajor        ("fmajor"        ,fmajor_p        ,2,2,2,nflav,ncol,nlay);
  umgInt4d  jeta          ("jeta"          ,jeta_p          ,2,    nflav,ncol,nlay);
  umgBool2d tropo         ("tropo"         ,tropo_p         ,ncol,nlay);
  umgInt2d  jtemp         ("jtemp"         ,jtemp_p         ,ncol,nlay);
  umgInt2d  jpress        ("jpress"        ,jpress_p        ,ncol,nlay);
  umgInt1d  gpoint_bands  ("gpoint_bands"  ,gpoint_bands_p  ,ngpt);
  umgInt2d  band_lims_gpt ("band_lims_gpt" ,band_lims_gpt_p ,2,nbnd);
  umgReal4d pfracin       ("pfracin"       ,pfracin_p       ,ngpt,neta,npres+1,ntemp);
  umgReal2d totplnk       ("totplnk"       ,totplnk_p       ,nPlanckTemp,nbnd);
  umgInt2d  gpoint_flavor ("gpoint_flavor" ,gpoint_flavor_p ,2,ngpt);
  umgReal2d sfc_src       ("sfc_src"       ,sfc_src_p       ,ngpt,     ncol);
  umgReal3d lay_src       ("lay_src"       ,lay_src_p       ,ngpt,nlay,ncol);
  umgReal3d lev_src_inc   ("lev_src_inc"   ,lev_src_inc_p   ,ngpt,nlay,ncol);
  umgReal3d lev_src_dec   ("lev_src_dec"   ,lev_src_dec_p   ,ngpt,nlay,ncol);

  real3d pfrac          ("pfrac"          ,ngpt,nlay,ncol);
  real3d planck_function("planck_function",nbnd,nlay+1,ncol);
  real1d one("one",2);
  one(1) = 1;
  one(2) = 1;

  // Calculation of fraction of band's Planck irradiance associated with each g-point
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  yakl::parallel_for_cpu_serial( Bounds<3>({1,ncol},{1,nlay},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int icol, ilay, igpt;
    yakl::storeIndices( indices , icol,ilay,igpt );

    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
    int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
    pfrac(igpt,ilay,icol) = 
      // interpolation in temperature, pressure, and eta
      interpolate3D(one, fmajor.slice(COLON,COLON,COLON,iflav,icol,ilay), pfracin, 
                    igpt, jeta.slice(COLON,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp);
  });

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  // for (int icol=1; icol<=ncol; icol++) {
  yakl::parallel_for_cpu_serial( Bounds<1>({1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int icol = indices[0];
    auto planck_function_slice = planck_function.slice(COLON,1,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });
  //
  // Map to g-points
  //
  // for (int igpt=1; igpt<=ngpt; igpt++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  yakl::parallel_for_cpu_serial( Bounds<2>({1,ngpt},{1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int igpt, icol;
    yakl::storeIndices( indices , igpt,icol );

    sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 1, icol);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  yakl::parallel_for_cpu_serial( Bounds<2>({1,ncol},{1,nlay}) , YAKL_LAMBDA ( int indices[] ) {
    int icol, ilay;
    yakl::storeIndices( indices , icol,ilay );

    // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
    auto planck_function_slice = planck_function.slice(COLON,ilay,icol); // Necessary to create a temporary because we're writing to it
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
  yakl::parallel_for_cpu_serial( Bounds<3>({1,ncol},{1,nlay},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int icol, ilay, igpt;
    yakl::storeIndices( indices , icol,ilay,igpt );

    lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
  });

  // compute level source irradiances for each g-point, one each for upward and downward paths
  // for (int icol=1; icol<=ncol; icol++) {
  yakl::parallel_for_cpu_serial( Bounds<1>({1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int icol = indices[0];
    auto planck_function_slice = planck_function.slice(COLON,1,icol); // Necessary to create a temporary because we're writing to it
    interpolate1D(tlev(icol,1), temp_ref_min, totplnk_delta, totplnk, planck_function_slice,nPlanckTemp,nbnd);
  });

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=2; ilay<=nlay+1; ilay++) {
  yakl::parallel_for_cpu_serial( Bounds<2>({1,ncol},{2,nlay+1}) , YAKL_LAMBDA ( int indices[] ) {
    int icol, ilay;
    yakl::storeIndices( indices , icol,ilay );

    auto planck_function_slice = planck_function.slice(COLON,ilay,icol); // Necessary to create a temporary because we're writing to it
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
  yakl::parallel_for_cpu_serial( Bounds<3>({1,ncol},{1,nlay},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int icol, ilay, igpt;
    yakl::storeIndices( indices , icol,ilay,igpt );

    lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
    lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
    if (icol < ncol) {
      lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
      lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
    }
  });

}



// ----------------------------------------------------------
//
// compute Rayleigh scattering optical depths
//
extern "C" void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta,
                                     int npres, int ntemp, int *gpoint_flavor_p, int *band_lims_gpt_p, 
                                     real *krayl_p, int idx_h2o, real *col_dry_p, real *col_gas_p,
                                     real *fminor_p, int *jeta_p, bool *tropo_p, int *jtemp_p,
                                     real *tau_rayleigh_p) {

  umgInt2d  gpoint_flavor("gpoint_flavor",gpoint_flavor_p,2,ngpt);
  umgInt2d  band_lims_gpt("band_lims_gpt",band_lims_gpt_p,2,nbnd);
  umgReal4d krayl        ("krayl"        ,krayl_p        ,ngpt,neta,ntemp,2);
  umgReal2d col_dry      ("col_dry"      ,col_dry_p      ,ncol,nlay);
  umgReal2d col_gas      ("col_gas"      ,col_gas_p      ,{1,ncol},{1,nlay},{0,ngas});
  umgReal5d fminor       ("fminor"       ,fminor_p       ,2,2,nflav,ncol,nlay);
  umgInt4d  jeta         ("jeta"         ,jeta_p         ,2,  nflav,ncol,nlay);
  umgBool2d tropo        ("tropo"        ,tropo_p        ,ncol,nlay);
  umgInt2d  jtemp        ("jtemp"        ,jtemp_p        ,ncol,nlay);
  umgReal3d tau_rayleigh ("tau_rayleigh" ,tau_rayleigh_p ,ngpt,nlay,ncol);

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  yakl::parallel_for_cpu_serial( Bounds<3>({1,nlay},{1,ncol},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int ilay, icol, igpt;
    yakl::storeIndices( indices , ilay,icol,igpt );

    int itropo = merge(1,2,tropo(icol,ilay)); // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int iflav = gpoint_flavor(itropo, igpt);
    real k = interpolate2D(fminor.slice(COLON,COLON,iflav,icol,ilay), 
                           krayl.slice(COLON,COLON,COLON,itropo),      
                           igpt, jeta.slice(COLON,iflav,icol,ilay), jtemp(icol,ilay), ngpt, neta, ntemp);
    tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay));
  });
}



// ----------------------------------------------------------
//
// compute minor species optical depths
//
void gas_optical_depths_minor(int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2d &gpt_flv,
                              real3d &kminor, int2d &minor_limits_gpt, bool1d &minor_scales_with_density,
                              bool1d &scale_by_complement, int1d &idx_minor, int1d &idx_minor_scaling,
                              int1d &kminor_start, real2d &play, real2d &tlay, real3d &col_gas, real5d &fminor,
                              int4d &jeta, int2d &layer_limits, int2d &jtemp, real3d &tau) {
  real constexpr PaTohPa = 0.01;

  int extent = nminor;

  // Find the largest number of g-points per band
  int max_gpt_diff = minor_limits_gpt(2,1) - minor_limits_gpt(1,1);
  for (int i=2; i<=nminor; i++) {
    max_gpt_diff = max( max_gpt_diff , minor_limits_gpt(2,i) - minor_limits_gpt(1,i) );
  }

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt0=0; igpt0<=max_gpt_diff; igpt0++) {
  yakl::parallel_for_cpu_serial( Bounds<3>({1,nlay},{1,ncol},{0,max_gpt_diff}) , YAKL_LAMBDA ( int indices[] ) {
    int ilay, icol, igpt0;
    yakl::storeIndices( indices , ilay,icol,igpt0 );
    
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
            real vmr_fact = 1._wp / mycol_gas_0;
            real dry_fact = 1._wp / (1._wp + mycol_gas_h2o * vmr_fact);
            // scale by density of special gas
            if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
              scaling = scaling * (1._wp - mycol_gas_imnr * vmr_fact * dry_fact);
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

          real tau_minor = 0._wp;
          int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
          int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
          real kminor_loc = interpolate2D(fminor.slice(COLON,COLON,iflav,icol,ilay), kminor, minor_loc,  
                                          jeta.slice(COLON,iflav,icol,ilay), myjtemp, nminork, neta, ntemp);
          tau_minor = kminor_loc * scaling;

          yakl::atomicAdd( tau(igpt,ilay,icol) , tau_minor );
        }  // igpt <= gptE
      }
    }
  });
}



// --------------------------------------------------------------------------------------
//
// compute minor species optical depths
//
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2d &gpoint_flavor, int2d &band_lims_gpt, real4d &kmajor,                         
                              real4d &col_mix, real6d &fmajor, int4d &jeta, bool2d &tropo, 
                              int2d &jtemp, int2d &jpress, real3d &tau) {
  // optical depth calculation for major species
  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     // optical depth calculation for major species
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  yakl::parallel_for_cpu_serial( Bounds<3>({1,nlay},{1,ncol},{1,ngpt}) , YAKL_LAMBDA ( int indices[] ) {
    int ilay, icol, igpt;
    yakl::storeIndices( indices , ilay,icol,igpt );

    // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
    int itropo = merge(1,2,tropo(icol,ilay));  // WS: moved inside innermost loop

    // binary species parameter (eta) and col_mix depend on band flavor
    int iflav = gpoint_flavor(itropo, igpt);
    real tau_major = 
      // interpolation in temperature, pressure, and eta
      interpolate3D(col_mix.slice(COLON,iflav,icol,ilay), 
                    fmajor.slice(COLON,COLON,COLON,iflav,icol,ilay), kmajor, 
                    igpt, jeta.slice(COLON,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp);
    tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_major;
  });
}



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
                                       int *jtemp_p, int *jpress_p, real *tau_p) {

  umgInt2d   gpoint_flavor                  ( "gpoint_flavor                  " , gpoint_flavor_p                   ,   2,ngpt                      );
  umgInt2d   band_lims_gpt                  ( "band_lims_gpt                  " , band_lims_gpt_p                   ,   2,nbnd                      );
  umgReal4d  kmajor                         ( "kmajor                         " , kmajor_p                          ,   ngpt,neta,npres+1,ntemp     );
  umgReal3d  kminor_lower                   ( "kminor_lower                   " , kminor_lower_p                    ,   nminorklower,neta,ntemp     );
  umgReal3d  kminor_upper                   ( "kminor_upper                   " , kminor_upper_p                    ,   nminorkupper,neta,ntemp     );
  umgInt2d   minor_limits_gpt_lower         ( "minor_limits_gpt_lower         " , minor_limits_gpt_lower_p          ,   2,nminorlower               );
  umgInt2d   minor_limits_gpt_upper         ( "minor_limits_gpt_upper         " , minor_limits_gpt_upper_p          ,   2,nminorupper               );
  umgBool1d  minor_scales_with_density_lower( "minor_scales_with_density_lower" , minor_scales_with_density_lower_p ,     nminorlower               );
  umgBool1d  minor_scales_with_density_upper( "minor_scales_with_density_upper" , minor_scales_with_density_upper_p ,     nminorupper               );
  umgBool1d  scale_by_complement_lower      ( "scale_by_complement_lower      " , scale_by_complement_lower_p       ,     nminorlower               );
  umgBool1d  scale_by_complement_upper      ( "scale_by_complement_upper      " , scale_by_complement_upper_p       ,     nminorupper               );
  umgInt1d   idx_minor_lower                ( "idx_minor_lower                " , idx_minor_lower_p                 ,     nminorlower               );
  umgInt1d   idx_minor_upper                ( "idx_minor_upper                " , idx_minor_upper_p                 ,     nminorupper               );
  umgInt1d   idx_minor_scaling_lower        ( "idx_minor_scaling_lower        " , idx_minor_scaling_lower_p         ,     nminorlower               );
  umgInt1d   idx_minor_scaling_upper        ( "idx_minor_scaling_upper        " , idx_minor_scaling_upper_p         ,     nminorupper               );
  umgInt1d   kminor_start_lower             ( "kminor_start_lower             " , kminor_start_lower_p              ,     nminorlower               );
  umgInt1d   kminor_start_upper             ( "kminor_start_upper             " , kminor_start_upper_p              ,     nminorupper               );
  umgBool2d  tropo                          ( "tropo                          " , tropo_p                           ,   ncol,nlay                   );
  umgReal4d  col_mix                        ( "col_mix                        " , col_mix_p                         ,   2,    nflav,ncol,nlay       );
  umgReal6d  fmajor                         ( "fmajor                         " , fmajor_p                          ,   2,2,2,nflav,ncol,nlay       );
  umgReal5d  fminor                         ( "fminor                         " , fminor_p                          ,   2,2,  nflav,ncol,nlay       );
  umgReal2d  play                           ( "play                           " , play_p                            ,               ncol,nlay       );
  umgReal2d  tlay                           ( "tlay                           " , tlay_p                            ,               ncol,nlay       );
  umgReal3d  col_gas                        ( "col_gas                        " , col_gas_p                         , {1,ncol},{1,nlay},{0,ngas}    );
  umgInt4d   jeta                           ( "jeta                           " , jeta_p                            ,   2,    nflav,ncol,nlay       );
  umgInt2d   jtemp                          ( "jtemp                          " , jtemp_p                           ,               ncol,nlay       );
  umgInt2d   jpress                         ( "jpress                         " , jpress_p                          ,               ncol,nlay       );
  umgReal3d  tau                            ( "tau                            " , tau_p                             ,   ngpt,nlay,ncol              );

  int2d itropo_lower("itropo_lower",ncol,2);
  int2d itropo_upper("itropo_upper",ncol,2);

  int huge  = std::numeric_limits<int>::max();
  int small = std::numeric_limits<int>::min();

  // ---------------------
  // Layer limits of upper, lower atmospheres
  // ---------------------
  bool top_at_1 = play(1,1) < play(1, nlay);

  if(top_at_1) {

    // for (int icol=1; icol<=ncol; icol++){
    yakl::parallel_for_cpu_serial( Bounds<1>({1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
      int icol = indices[0];

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
    yakl::parallel_for_cpu_serial( Bounds<1>({1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
      int icol = indices[0];

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
  gas_optical_depths_major(   
   ncol,nlay,nbnd,ngpt,       
   nflav,neta,npres,ntemp,    
   gpoint_flavor,             
   band_lims_gpt,             
   kmajor,                    
   col_mix,fmajor,            
   jeta,tropo,jtemp,jpress,   
   tau);
  // ---------------------
  // Minor Species - lower
  // ---------------------
  int idx_tropo = 1;
  gas_optical_depths_minor(     
    ncol,nlay,ngpt,             
    ngas,nflav,ntemp,neta,      
    nminorlower,nminorklower,   
    idx_h2o,idx_tropo,          
    gpoint_flavor,              
    kminor_lower,               
    minor_limits_gpt_lower,     
    minor_scales_with_density_lower,
    scale_by_complement_lower,  
    idx_minor_lower,            
    idx_minor_scaling_lower,    
    kminor_start_lower,         
    play, tlay,                 
    col_gas,fminor,jeta,        
    itropo_lower,jtemp,         
    tau);
  // ---------------------
  // Minor Species - upper
  // ---------------------
  idx_tropo = 2;
  gas_optical_depths_minor(     
    ncol,nlay,ngpt,             
    ngas,nflav,ntemp,neta,      
    nminorupper,nminorkupper,   
    idx_h2o,idx_tropo,          
    gpoint_flavor,              
    kminor_upper,               
    minor_limits_gpt_upper,     
    minor_scales_with_density_upper, 
    scale_by_complement_upper,  
    idx_minor_upper,            
    idx_minor_scaling_upper,    
    kminor_start_upper,         
    play, tlay,                 
    col_gas,fminor,jeta,        
    itropo_upper,jtemp,         
    tau);

}





