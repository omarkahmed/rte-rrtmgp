
#include "mo_gas_optics_kernels.h"
#include <limits>

using yakl::dim2;

extern "C" void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp,
                              int *flavor_p, real *press_ref_log_p, real *temp_ref_p, real press_ref_log_delta,
                              real temp_ref_min, real temp_ref_delta, real press_ref_trop_log,
                              real *vmr_ref_p, real *play_p, real *tlay_p, real *col_gas_p, int *jtemp_p,
                              real *fmajor_p, real *fminor_p, real *col_mix_p, bool *tropo_p, int *jeta_p,
                              int *jpress_p) {

  umgInt2d  flavor       ("flavor"       ,flavor_p       ,2,nflav);
  umgReal1d press_ref_log("press_ref_log",press_ref_log_p,npres);
  umgReal1d temp_ref     ("temp_ref"     ,temp_ref_p     ,ntemp);
  umgReal2d vmr_ref      ("vmr_ref"      ,vmr_ref_p      , dim2(1,2) , dim2(0,ngas) , dim2(1,ntemp) );
  umgReal2d play         ("play"         ,play_p         ,ncol,nlay);
  umgReal2d tlay         ("tlay"         ,tlay_p         ,ncol,nlay);
  umgReal3d col_gas      ("col_gas"      ,col_gas_p      , dim2(1,ncol) , dim2(1,nlay) , dim2(0,ngas) );
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
          yakl::SArray<int,2> igases;

          // itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
          int itropo = merge(1,2,tropo(icol,ilay));
          igases(0) = flavor(1,iflav);
          igases(1) = flavor(2,iflav);

          // compute interpolation fractions needed for lower, then upper reference temperature level
          // compute binary species parameter (eta) for flavor and temperature and
          //  associated interpolation index and factors
          real ratio_eta_half = vmr_ref(itropo,igases(0),(jtemp(icol,ilay)+itemp-1)) / 
                                vmr_ref(itropo,igases(1),(jtemp(icol,ilay)+itemp-1));
          col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases(0)) + ratio_eta_half * col_gas(icol,ilay,igases(1));
          real eta = merge(col_gas(icol,ilay,igases(0)) / col_mix(itemp,iflav,icol,ilay), 0.5_wp, 
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

