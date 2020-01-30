
#include "mo_gas_optics_kernels.h"
#include <limits>

extern "C" void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp,
                              int *flavor_p, real *press_ref_log_p, real *temp_ref_p, real press_ref_log_delta,
                              real temp_ref_min, real temp_ref_delta, real press_ref_trop_log,
                              real *vmr_ref_p, real *play_p, real *tlay_p, real *col_gas_p, int *jtemp_p,
                              real *fmajor_p, real *fminor_p, real *col_mix_p, bool *tropo_p, int *jeta_p,
                              int *jpress_p) {
  umgInt2d  flavor       ("flavor"       ,flavor_p       ,nflav,2);
  umgReal1d press_ref_log("press_ref_log",press_ref_log_p,npres);
  umgReal1d temp_ref     ("temp_ref"     ,temp_ref_p     ,ntemp);
  umgReal2d vmr_ref      ("vmr_ref"      ,vmr_ref_p      ,ntemp,ngas+1,2);
  umgReal2d play         ("play"         ,play_p         ,nlay,ncol);
  umgReal2d tlay         ("tlay"         ,tlay_p         ,nlay,ncol);
  umgReal3d col_gas      ("col_gas"      ,col_gas_p      ,ngas+1,nlay,ncol);
  umgInt2d  jtemp        ("jtemp"        ,jtemp_p        ,nlay,ncol);
  umgInt2d  jpress       ("jpress"       ,jpress_p       ,nlay,ncol);
  umgBool2d tropo        ("tropo"        ,tropo_p        ,nlay,ncol);
  umgInt4d  jeta         ("jeta"         ,jeta_p         ,nlay,ncol,nflav,2);
  umgReal4d col_mix      ("col_mix"      ,col_mix_p      ,nlay,ncol,nflav,2);
  umgReal6d fmajor       ("fmajor"       ,fmajor_p       ,nlay,ncol,nflav,2,2,2);
  umgReal5d fminor       ("fminor"       ,fminor_p       ,nlay,ncol,nflav,2,2);

  int const off_gas = 1; // For vmr_ref and col_gas

  real2d ftemp ( "ftemp"  , nlay , ncol );
  real2d fpress( "fpress" , nlay , ncol );

  real tiny = std::numeric_limits<real>::min();

  for (int ilay=0; ilay<nlay; ilay++) {
    for (int icol=0; icol<ncol; ncol++) {
      // index and factor for temperature interpolation
      jtemp(ilay,icol) = (int) ((tlay(ilay,icol) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
      jtemp(ilay,icol) = min(ntemp - 1, max(1, jtemp(ilay,icol)));                                      // limit the index range
      ftemp(ilay,icol) = (tlay(ilay,icol) - temp_ref(jtemp(ilay,icol)-1)) / temp_ref_delta;

      // index and factor for pressure interpolation
      float locpress = 1. + (log(play(ilay,icol)) - press_ref_log(1)) / press_ref_log_delta;
      jpress(ilay,icol) = min(npres-1, max(1, (int) locpress));
      fpress(ilay,icol) = locpress - (real) jpress(ilay,icol);

      // determine if in lower or upper part of atmosphere
      tropo(ilay,icol) = log(play(ilay,icol)) > press_ref_trop_log;
    }
  }

  for (int ilay=0; ilay<nlay; ilay++) {
    for (int icol=0; icol<ncol; icol++) {
      for (int iflav=0; iflav<nflav; iflav++) {  // loop over implemented combinations of major species
        for (int itemp=1; itemp<=ntemp; itemp++) {
          yakl::SArray<int,2> igases;

          // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
          int itropo = tropo(ilay,icol) ? 0 : 1;
          igases(0) = flavor(iflav,0);
          igases(1) = flavor(iflav,1);

          // compute interpolation fractions needed for lower, then upper reference temperature level
          // compute binary species parameter (eta) for flavor and temperature and
          //  associated interpolation index and factors
          real vmr_ref1 = vmr_ref((jtemp(icol,ilay)+itemp-2),off_gas+igases(0),itropo);
          real vmr_ref2 = vmr_ref((jtemp(icol,ilay)+itemp-2),off_gas+igases(1),itropo);
          real ratio_eta_half = vmr_ref1 / vmr_ref2;
          col_mix(ilay,icol,iflav,itemp-1) = col_gas(off_gas+igases(0),ilay,icol) + ratio_eta_half * col_gas(off_gas+igases(1),ilay,icol);
          real eta = col_mix(ilay,icol,iflav,itemp-1) > 2. * tiny ? col_gas(off_gas+igases(1),ilay,icol) / col_mix(ilay,icol,iflav,itemp-1) : 0.5;
          real loceta = eta * (real) (neta-1);
          jeta(ilay,icol,iflav,itemp-1) = min((int)loceta+1, neta-1);
          real feta = fmod( loceta , 1.0 );
          // compute interpolation fractions needed for minor species
          real ftemp_term = ( (2.0-itemp) + (2.0*itemp-3.0) * ftemp(ilay,icol));
          fminor(ilay,icol,iflav,itemp-1,0) = (1.-feta) * ftemp_term;
          fminor(ilay,icol,iflav,itemp-1,1) =     feta  * ftemp_term;
          // compute interpolation fractions needed for major species
          fmajor(ilay,icol,iflav,itemp-1,0,0) = (1.-fpress(ilay,icol)) * fminor(ilay,icol,iflav,itemp-1,0);
          fmajor(ilay,icol,iflav,itemp-1,0,1) = (1.-fpress(ilay,icol)) * fminor(ilay,icol,iflav,itemp-1,1);
          fmajor(ilay,icol,iflav,itemp-1,1,0) =     fpress(ilay,icol)  * fminor(ilay,icol,iflav,itemp-1,0);
          fmajor(ilay,icol,iflav,itemp-1,1,1) =     fpress(ilay,icol)  * fminor(ilay,icol,iflav,itemp-1,1);
        }
      }
    }
  }

}
