
#include "mo_optical_props_kernels.h"



// ---------------------------------
// increment 2stream by 2stream
extern "C" void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                             real *tau1_p, real *ssa1_p, real *g1_p,
                                             real *tau2_p, real *ssa2_p, real *g2_p,
                                             int nbnd, int *gpt_lims_p) {
  umgReal3d tau1     ("tau1"     ,tau1_p     ,ncol,nlay,ngpt);
  umgReal3d tau2     ("tau2"     ,tau2_p     ,ncol,nlay,nbnd);
  umgReal3d ssa1     ("ssa1"     ,ssa1_p     ,ncol,nlay,ngpt);
  umgReal3d ssa2     ("ssa2"     ,ssa2_p     ,ncol,nlay,nbnd);
  umgReal3d g1       ("g1"       ,g1_p       ,ncol,nlay,ngpt);
  umgReal3d g2       ("g2"       ,g2_p       ,ncol,nlay,nbnd);
  umgInt2d  gpt_lims ("gpt_lims" ,gpt_lims_p ,2,nbnd);

  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1 , ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay},{1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );

    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1,ibnd) && igpt <= gpt_lims(2,ibnd) ) {
        // t=tau1 + tau2
        real tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
        // w=(tau1*ssa1 + tau2*ssa2) / t
        real tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) + 
                         tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
        g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(icol,ilay,igpt) + 
                              tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * g2(icol,ilay,ibnd)) / max(eps,tauscat12);
        ssa1(icol,ilay,igpt) = tauscat12 / max(eps,tau12);
        tau1(icol,ilay,igpt) = tau12;
      }
    }
  });
}



// ---------------------------------
//
// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
//
// ---------------------------------
extern "C" void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real *tau1_p, real *tau2_p,             
                                             int nbnd, int *gpt_lims_p) {
  umgReal3d tau1    ("tau1"    ,tau1_p    ,ncol,nlay,ngpt);
  umgReal3d tau2    ("tau2"    ,tau2_p    ,ncol,nlay,nbnd);
  umgInt2d  gpt_lims("gpt_lims",gpt_lims_p,2,nbnd);

  // do igpt = 1 , ngpt
  //   do ilay = 1 , nlay
  //     do icol = 1 , ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay},{1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );

    for (int ibnd=1; ibnd<=nbnd; ibnd++) {
      if (igpt >= gpt_lims(1,ibnd) && igpt <= gpt_lims(2,ibnd) ) {
        tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
      }
    }
  });
}



// ---------------------------------
// Delta-scale
//   f = g*g
//
extern "C" void delta_scale_2str_k(int ncol, int nlay, int ngpt, real *tau_p, real *ssa_p, real *g_p) {
  umgReal3d tau("tau",tau_p,ncol,nlay,ngpt);
  umgReal3d ssa("ssa",ssa_p,ncol,nlay,ngpt);
  umgReal3d g  ("g"  ,g_p  ,ncol,nlay,ngpt);

  real eps = 3*std::numeric_limits<real>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({1,ngpt},{1,nlay},{1,ncol}) , YAKL_LAMBDA ( int indices[] ) {
    int igpt, ilay, icol;
    storeIndices( indices , igpt,ilay,icol );

    if (tau(icol,ilay,igpt) > eps) {
      real f  = g  (icol,ilay,igpt) * g  (icol,ilay,igpt);
      real wf = ssa(icol,ilay,igpt) * f;
      tau(icol,ilay,igpt) = (1._wp - wf) * tau(icol,ilay,igpt);
      ssa(icol,ilay,igpt) = (ssa(icol,ilay,igpt) - wf) / (1.0_wp - wf);
      g  (icol,ilay,igpt) = (g  (icol,ilay,igpt) -  f) / (1.0_wp -  f);
    }
  });
}




