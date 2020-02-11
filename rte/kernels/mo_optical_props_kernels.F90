! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Description:  Addition of optical properties -- the first set are incremented by the second set.
!   There are three possible representations of optical properties (scalar = optical depth only;
!   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
!   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
!   each choice of representation on the left hand side times three representations of the
!   optical properties to be added.
!   There are two sets of these nine routines. In the first the two sets of optical
!   properties are defined at the same spectral resolution. There is also a set of routines
!   to add properties defined at lower spectral resolution to a set defined at higher spectral
!   resolution (adding properties defined by band to those defined by g-point)

module mo_optical_props_kernels
  use, intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  implicit none
  public


  interface
    subroutine inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt, &
                                                 tau1, ssa1, g1,   &
                                                 tau2, ssa2, g2,   &
                                                 nbnd, gpt_lims) bind(C, name="inc_2stream_by_2stream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2, g2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine

    subroutine inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                            tau1,             &
                                            tau2,             &
                                            nbnd, gpt_lims) bind(C, name="inc_1scalar_by_1scalar_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine

    subroutine delta_scale_2str_f_k(ncol, nlay, ngpt, tau, ssa, g, f) &
        bind(C, name="delta_scale_2str_f_k")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g
      real(wp), dimension(ncol, nlay, ngpt), intent(in   ) ::  f
    end subroutine

    subroutine increment_1scalar_by_1scalar(ncol, nlay, ngpt, &
                                                 tau1,             &
                                                 tau2) bind(C, name="increment_1scalar_by_1scalar")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in  ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2
    end subroutine

    subroutine increment_1scalar_by_2stream(ncol, nlay, ngpt, &
                                                 tau1,             &
                                                 tau2, ssa2) bind(C, name="increment_1scalar_by_2stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
    end subroutine increment_1scalar_by_2stream

    subroutine increment_1scalar_by_nstream(ncol, nlay, ngpt, &
                                                 tau1,             &
                                                 tau2, ssa2) bind(C, name="increment_1scalar_by_nstream")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
    end subroutine increment_1scalar_by_nstream

    subroutine increment_2stream_by_1scalar(ncol, nlay, ngpt, &
                                                 tau1, ssa1,       &
                                                 tau2) bind(C, name="increment_2stream_by_1scalar")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2
    end subroutine increment_2stream_by_1scalar

    subroutine increment_2stream_by_2stream(ncol, nlay, ngpt, &
                                                 tau1, ssa1, g1,   &
                                                 tau2, ssa2, g2) bind(C, name="increment_2stream_by_2stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2, g2
    end subroutine increment_2stream_by_2stream

    subroutine increment_2stream_by_nstream(ncol, nlay, ngpt, nmom2, &
                                                 tau1, ssa1, g1,          &
                                                 tau2, ssa2, p2) bind(C, name="increment_2stream_by_nstream")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt, nmom2
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
      real(wp), dimension(nmom2, &
                          ncol,nlay,ngpt), intent(in   ) :: p2
    end subroutine increment_2stream_by_nstream

    subroutine increment_nstream_by_1scalar(ncol, nlay, ngpt, &
                                                 tau1, ssa1,       &
                                                 tau2) bind(C, name="increment_nstream_by_1scalar")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2
    end subroutine increment_nstream_by_1scalar

    subroutine increment_nstream_by_2stream(ncol, nlay, ngpt, nmom1, &
                                                 tau1, ssa1, p1,          &
                                                 tau2, ssa2, g2) bind(C, name="increment_nstream_by_2stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                       intent(in   ) :: ncol, nlay, ngpt, nmom1
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(nmom1, &
                          ncol,nlay,ngpt), intent(inout) :: p1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2, g2
    end subroutine increment_nstream_by_2stream

    subroutine increment_nstream_by_nstream(ncol, nlay, ngpt, nmom1, nmom2, &
                                                 tau1, ssa1, p1,                 &
                                                 tau2, ssa2, p2) bind(C, name="increment_nstream_by_nstream")
      use mo_rte_kind, only: wp, wl
      integer,  value,                      intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(nmom1, &
                          ncol,nlay,ngpt), intent(inout) :: p1
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
      real(wp), dimension(nmom2, &
                          ncol,nlay,ngpt), intent(in   ) :: p2
    end subroutine increment_nstream_by_nstream

    subroutine inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, &
                                                 tau1,             &
                                                 tau2, ssa2,       &
                                                 nbnd, gpt_lims) bind(C, name="inc_1scalar_by_2stream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_1scalar_by_2stream_bybnd

    subroutine inc_1scalar_by_nstream_bybnd(ncol, nlay, ngpt, &
                                                 tau1,             &
                                                 tau2, ssa2,       &
                                                 nbnd, gpt_lims) bind(C, name="inc_1scalar_by_nstream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_1scalar_by_nstream_bybnd

    subroutine inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                                 tau1, ssa1,       &
                                                 tau2,             &
                                                 nbnd, gpt_lims) bind(C, name="inc_2stream_by_1scalar_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_2stream_by_1scalar_bybnd

    subroutine inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, nmom2, &
                                                 tau1, ssa1, g1,          &
                                                 tau2, ssa2, p2,          &
                                                 nbnd, gpt_lims) bind(C, name="inc_2stream_by_nstream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nmom2, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
      real(wp), dimension(nmom2, &
                          ncol,nlay,nbnd), intent(in   ) :: p2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_2stream_by_nstream_bybnd
  
    subroutine inc_nstream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                                 tau1, ssa1,       &
                                                 tau2,             &
                                                 nbnd, gpt_lims) bind(C, name="inc_nstream_by_1scalar_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_nstream_by_1scalar_bybnd

    subroutine inc_nstream_by_2stream_bybnd(ncol, nlay, ngpt, nmom1, &
                                                 tau1, ssa1, p1,          &
                                                 tau2, ssa2, g2,          &
                                                 nbnd, gpt_lims) bind(C, name="inc_nstream_by_2stream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nmom1, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(nmom1, &
                          ncol,nlay,ngpt), intent(inout) :: p1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2, g2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_nstream_by_2stream_bybnd

    subroutine inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, nmom1, nmom2, &
                                                 tau1, ssa1, p1,                 &
                                                 tau2, ssa2, p2,                 &
                                                 nbnd, gpt_lims) bind(C, name="inc_nstream_by_nstream_bybnd")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2, nbnd
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
      real(wp), dimension(nmom1, &
                          ncol,nlay,ngpt), intent(inout) :: p1
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
      real(wp), dimension(nmom2, &
                          ncol,nlay,nbnd), intent(in   ) :: p2
      integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims ! Starting and ending gpoint for each band
    end subroutine inc_nstream_by_nstream_bybnd

    subroutine extract_subset_dim1_3d(ncol, nlay, ngpt, array_in, colS, colE, array_out) &
      bind (C, name="extract_subset_dim1_3d")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(in ) :: array_in
      integer, value,                      intent(in ) :: colS, colE
      real(wp), dimension(colE-colS+1,&
                               nlay,ngpt), intent(out) :: array_out
    end subroutine extract_subset_dim1_3d

    subroutine extract_subset_dim2_4d(nmom, ncol, nlay, ngpt, array_in, colS, colE, array_out) &
      bind (C, name="extract_subset_dim2_4d")
      use mo_rte_kind, only: wp, wl
      integer, value,                           intent(in ) :: nmom, ncol, nlay, ngpt
      real(wp), dimension(nmom,ncol,nlay,ngpt), intent(in ) :: array_in
      integer, value,                           intent(in ) :: colS, colE
      real(wp), dimension(nmom,colE-colS+1,&
                                    nlay,ngpt), intent(out) :: array_out
    end subroutine extract_subset_dim2_4d

    subroutine extract_subset_absorption_tau(ncol, nlay, ngpt, tau_in, ssa_in, &
                                                  colS, colE, tau_out)              &
      bind (C, name="extract_subset_absorption_tau")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(in ) :: tau_in, ssa_in
      integer,                             intent(in ) :: colS, colE
      real(wp), dimension(colE-colS+1,&
                               nlay,ngpt), intent(out) :: tau_out
    end subroutine extract_subset_absorption_tau

    subroutine delta_scale_2str_k(ncol, nlay, ngpt, tau, ssa, g) bind(C, name="delta_scale_2str_k")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt
      real(wp), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g
    end subroutine

  end interface

  interface delta_scale_2str_kernel
    procedure delta_scale_2str_f_k, delta_scale_2str_k
  end interface

  interface extract_subset
    procedure extract_subset_dim1_3d, extract_subset_dim2_4d
    procedure extract_subset_absorption_tau
  end interface extract_subset

end module mo_optical_props_kernels
