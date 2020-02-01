
#include "mo_gas_optics_kernels.h"
#include <limits>

using yakl::bnd;
using yakl::COLON;



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

  for (int ilay=1; ilay<=nlay; ilay++) {
    for (int icol=1; icol<=ncol; icol++) {
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
    }
  }

  for (int ilay=1; ilay<=nlay; ilay++) {
    for (int icol=1; icol<=ncol; icol++) {
      for (int iflav=1; iflav<=nflav; iflav++) {   // loop over implemented combinations of major species
        for (int itemp=1; itemp<=2; itemp++) {
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
        }
      }
    }
  }
}



extern "C" void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real *tau_abs_p, real *tau_rayleigh_p,
                                         real *tau_p,  real *ssa_p, real *g_p) {

  umgReal3d tau_abs     ("tau_abs"     ,tau_abs_p     ,ngpt,nlay,ncol);
  umgReal3d tau_rayleigh("tau_rayleigh",tau_rayleigh_p,ngpt,nlay,ncol);
  umgReal3d tau         ("tau"         ,tau_p         ,ncol,nlay,ngpt);
  umgReal3d ssa         ("ssa"         ,ssa_p         ,ncol,nlay,ngpt);
  umgReal3d g           ("g"           ,g_p           ,ncol,nlay,ngpt);

  real tiny = std::numeric_limits<real>::min();

  for (int icol=1; icol<=ncol; icol++) {
    for (int ilay=1; ilay<=nlay; ilay++) {
      for (int igpt=1; igpt<=ngpt; igpt++) {
         real t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
         tau(icol,ilay,igpt) = t;
         g  (icol,ilay,igpt) = 0._wp;
         if(t > 2._wp * tiny) {
           ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
         } else {
           ssa(icol,ilay,igpt) = 0._wp;
         }
      }
    }
  }
}



// ----------------------------------------------------------
// interpolation in temperature, pressure, and eta
extern "C" real interpolate3D(real *scaling_p, real *fmajor_p, real *k_p, int igpt, int *jeta_p,
                              int jtemp, int jpress, int ngpt, int neta, int npres, int ntemp) {
  
  umgReal1d scaling("scaling",scaling_p,2);
  umgReal3d fmajor ("fmajor" ,fmajor_p ,2,2,2);
  umgReal4d k      ("k"      ,k_p      ,ngpt,neta,npres+1,ntemp);
  umgInt1d  jeta   ("jeta"   ,jeta_p   ,2);

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
real interpolate3D(real1d scaling, real3d fmajor, real4d k, int igpt, int1d jeta,
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
extern "C" real interpolate2D(real *fminor_p, real *k_p, int igpt, int *jeta_p, int jtemp,
                              int ngpt, int neta, int ntemp) {

  umgReal2d fminor("fminor",fminor_p,2,2);
  umgReal3d k     ("k"     ,k_p     ,ngpt,neta,ntemp);
  umgInt1d  jeta  ("jeta"  ,jeta_p  ,2);

  return fminor(1,1) * k(igpt, jeta(1)  , jtemp  ) + 
         fminor(2,1) * k(igpt, jeta(1)+1, jtemp  ) + 
         fminor(1,2) * k(igpt, jeta(2)  , jtemp+1) + 
         fminor(2,2) * k(igpt, jeta(2)+1, jtemp+1);
}
real interpolate2D(real2d fminor, real3d k, int igpt, int1d jeta, int jtemp,
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
extern "C" void interpolate1D(real val, real offset, real delta, real *table_p,
                              real *res_p, int tab_d1, int tab_d2) {

  umgReal2d table("table",table_p,tab_d1,tab_d2);
  umgReal1d res  ("res"  ,res_p  ,tab_d2);

  real val0 = (val - offset) / delta;
  real frac = val0 - int(val0); // get fractional part
  int index = min(tab_d1-1, max(1, (int)(val0)+1)); // limit the index range
  for (int i=1; i<=tab_d2; i++) {
    res(i) = table(index,i) + frac * (table(index+1,i) - table(index,i));
  }
}
void interpolate1D(real val, real offset, real delta, real2d table,
                   real2d res, int tab_d1, int tab_d2) {
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

  real3d pfrac          ("pfrac",ngpt,nlay,ncol);
  real3d planck_function("pfrac",nbnd,nlay+1,ncol);
  real1d one("one",2);
  one(1) = 1;
  one(2) = 1;

  // Calculation of fraction of band's Planck irradiance associated with each g-point
  for (int icol=1; icol<=ncol; icol++) {
    for (int ilay=1; ilay<=nlay; ilay++) {
      for (int igpt=1; igpt<=ngpt; igpt++) {
        // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        int itropo = merge(1,2,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
        int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
        pfrac(igpt,ilay,icol) = 
          // interpolation in temperature, pressure, and eta
          interpolate3D(one, fmajor.slice(COLON,COLON,COLON,iflav,icol,ilay), pfracin, 
                        igpt, jeta.slice(COLON,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp);
      }
    }
  }

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  for (int icol=1; icol<=ncol; icol++) {
    interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, planck_function.slice(COLON,1,icol),nPlanckTemp,nbnd);
  }
  //
  // Map to g-points
  //
  for (int igpt=1; igpt<=ngpt; igpt++) {
    for (int icol=1; icol<=ncol; icol++) {
      sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 1, icol);
    }
  }

  for (int icol=1; icol<=ncol; icol++) {
    for (int ilay=1; ilay<=nlay; ilay++) {
      // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
      interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function.slice(COLON,ilay,icol),nPlanckTemp,nbnd);
    }
  }
  //
  // Map to g-points
  //
  // Explicitly unroll a time-consuming loop here to increase instruction-level parallelism on a GPU
  // Helps to achieve higher bandwidth
  //
  for (int icol=1; icol<=ncol; icol+=2) {
    for (int ilay=1; ilay<=nlay; ilay++) {
      for (int igpt=1; igpt<=ngpt; igpt++) {
        lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
        if (icol < ncol) {
          lay_src(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,icol+1);
        }
      }
    }
  }

  // compute level source irradiances for each g-point, one each for upward and downward paths
  for (int icol=1; icol<=ncol; icol++) {
    interpolate1D(tlev(icol,1), temp_ref_min, totplnk_delta, totplnk, planck_function.slice(COLON,1,icol),nPlanckTemp,nbnd);
  }

  for (int icol=1; icol<=ncol; icol++) {
    for (int ilay=2; ilay<=nlay+1; ilay++) {
      interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk, planck_function.slice(COLON,ilay,icol),nPlanckTemp,nbnd);
    }
  }

  //
  // Map to g-points
  //
  // Same unrolling as mentioned before
  //
  for (int icol=1; icol<=ncol; icol+=2) {
    for (int ilay=1; ilay<=nlay; ilay++) {
      for (int igpt=1; igpt<=ngpt; igpt++) {
        lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
        lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
        if (icol < ncol) {
          lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
          lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
        }
      }
    }
  }
  // do icol = 1, ncol, 2
  //   do ilay = 1, nlay
  //     do igpt = 1, ngpt
  //       lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  )
  //       lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  )
  //       if (icol < ncol) then
  //       lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1)
  //       lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1)
  //       end if
  //     end do
  //   end do ! ilay
  // end do ! icol

}

