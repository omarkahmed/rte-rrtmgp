! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Numeric calculations for radiative transfer solvers.
!   Emission/absorption (no-scattering) calculations
!     solver for multi-angle Gaussian quadrature
!     solver for a single angle, calling
!       source function computation (linear-in-tau)
!       transport
!   Extinction-only calculation (direct solar beam)
!   Two-stream calculations
!     solvers for LW and SW with different boundary conditions and source functions
!       source function calculation for LW, SW
!       two-stream calculations for LW, SW (using different assumtions about phase function)
!       transport (adding)
!   Application of boundary conditions
!
! -------------------------------------------------------------------------------------------------
module mo_rte_solver_kernels
  use,  intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  implicit none
  private

  interface

    subroutine lw_source_noscat_stencil(ncol, nlay, ngpt, icol, ilay, igpt,                   &
                                        lay_source, lev_source_up, lev_source_dn, tau, trans, &
                                        source_dn, source_up) bind(C, name="lw_source_noscat_stencil")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in)   :: ncol, nlay, ngpt
      integer, value,                        intent(in)   :: icol, ilay, igpt ! Working point coordinates
      real(wp), dimension(ncol, nlay, ngpt), intent(in)   :: lay_source,    & ! Planck source at layer center
                                                             lev_source_up, & ! Planck source at levels (layer edges),
                                                             lev_source_dn, & !   increasing/decreasing layer index
                                                             tau,           & ! Optical path (tau/mu)
                                                             trans            ! Transmissivity (exp(-tau))
      real(wp), dimension(ncol, nlay, ngpt), intent(inout):: source_dn, source_up
    end subroutine

    subroutine apply_BC_0(ncol, nlay, ngpt, top_at_1, flux_dn) bind(C, name="apply_BC_0")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_dn          ! Flux to be used as input to solvers below
    end subroutine

    subroutine apply_BC_factor(ncol, nlay, ngpt, top_at_1, inc_flux, factor, flux_dn) bind(C, name="apply_BC_factor")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: inc_flux         ! Flux at top of domain
      real(wp), dimension(ncol            ), intent(in   ) :: factor           ! Factor to multiply incoming flux
      real(wp), dimension(ncol,nlay+1,ngpt), intent(out  ) :: flux_dn          ! Flux to be used as input to solvers below
    end subroutine

    subroutine apply_BC_gpt(ncol, nlay, ngpt, top_at_1, inc_flux, flux_dn) bind(C, name="apply_BC_gpt")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: inc_flux         ! Flux at top of domain
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_dn          ! Flux to be used as input to solvers below
    end subroutine

    subroutine adding(ncol, nlay, ngpt, top_at_1, &
                      albedo_sfc,           &
                      rdif, tdif,           &
                      src_dn, src_up, src_sfc, &
                      flux_up, flux_dn) bind(C, name="adding")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol       ,ngpt), intent(in   ) :: albedo_sfc
      real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: rdif, tdif
      real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: src_dn, src_up
      real(wp), dimension(ncol       ,ngpt), intent(in   ) :: src_sfc
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn
    end subroutine

    subroutine sw_source_2str(ncol, nlay, ngpt, top_at_1, Rdir, Tdir, Tnoscat, sfc_albedo, &
                              source_up, source_dn, source_sfc, flux_dn_dir) bind(C, name="sw_source_2str")
      use mo_rte_kind, only: wp, wl
      integer, value,                          intent(in   ) :: ncol, nlay, ngpt
      logical(wl), value,                      intent(in   ) :: top_at_1
      real(wp), dimension(ncol, nlay  , ngpt), intent(in   ) :: Rdir, Tdir, Tnoscat ! Layer reflectance, transmittance for diffuse radiation
      real(wp), dimension(ncol        , ngpt), intent(in   ) :: sfc_albedo          ! surface albedo for direct radiation
      real(wp), dimension(ncol, nlay  , ngpt), intent(  out) :: source_dn, source_up
      real(wp), dimension(ncol        , ngpt), intent(  out) :: source_sfc          ! Source function for upward radation at surface
      real(wp), dimension(ncol, nlay+1, ngpt), intent(inout) :: flux_dn_dir ! Direct beam flux
    end subroutine

    subroutine sw_two_stream(ncol, nlay, ngpt, mu0, tau, w0, g, &
                                  Rdif, Tdif, Rdir, Tdir, Tnoscat) bind (C, name="sw_two_stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in)  :: ncol, nlay, ngpt
      real(wp), dimension(ncol),           intent(in)  :: mu0
      real(wp), dimension(ncol,nlay,ngpt), intent(in)  :: tau, w0, g
      real(wp), dimension(ncol,nlay,ngpt), intent(out) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
    end subroutine

    subroutine lw_transport_noscat(ncol, nlay, ngpt, top_at_1, &
                                   tau, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                   radn_up, radn_dn) bind(C, name="lw_transport_noscat")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1   !
      real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                              trans      ! transmissivity = exp(-tau)
      real(wp), dimension(ncol       ,ngpt), intent(in   ) :: sfc_albedo ! Surface albedo
      real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: source_dn, &
                                                              source_up  ! Diffuse radiation emitted by the layer
      real(wp), dimension(ncol       ,ngpt), intent(in   ) :: source_sfc ! Surface source function [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn ! Radiances [W/m2-str]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up ! Radiances [W/m2-str]
    end subroutine

    subroutine sw_solver_2stream (ncol, nlay, ngpt, top_at_1, &
                                  tau, ssa, g, mu0,           &
                                  sfc_alb_dir, sfc_alb_dif,   &
                                  flux_up, flux_dn, flux_dir) bind (C, name="sw_solver_2stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau, &  ! Optical thickness,
                                                              ssa, &  ! single-scattering albedo,
                                                              g       ! asymmetry parameter []
      real(wp), dimension(ncol            ), intent(in   ) :: mu0     ! cosine of solar zenith angle
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_alb_dir, sfc_alb_dif
                                                                    ! Spectral albedo of surface to direct and diffuse radiation
      real(wp), dimension(ncol,nlay+1,ngpt), &
                                              intent(  out) :: flux_up ! Fluxes [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), &                        ! Downward fluxes contain boundary conditions
                                              intent(inout) :: flux_dn, flux_dir
    end subroutine

    subroutine lw_solver_noscat(ncol, nlay, ngpt, top_at_1, D, weight,                             &
                                tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                                radn_up, radn_dn) bind(C, name="lw_solver_noscat")
      use mo_rte_kind, only: wp, wl
      integer, value,             intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,         intent( in) :: top_at_1
      real(wp), dimension(ncol,       ngpt), intent( in) :: D            ! secant of propagation angle  []
      real(wp), value,                       intent( in) :: weight       ! quadrature weight
      real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
      real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
      ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
      ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
      ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
      real(wp), dimension(ncol,nlay,  ngpt), target, &
                                             intent( in) :: lev_source_inc, lev_source_dec
      real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
      real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn ! Radiances [W/m2-str]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up ! Radiances [W/m2-str]
    end subroutine

    subroutine lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                     tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_up, flux_dn) &
                                     bind (C, name="lw_solver_noscat_GaussQuad")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      integer, value,                        intent(in   ) :: nmus          ! number of quadrature angles
      real(wp), dimension(nmus),             intent(in   ) :: Ds, weights  ! quadrature secants, weights
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lev_source_inc
                                          ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                          ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lev_source_dec
                                                 ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis         ! Surface emissivity      []
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src          ! Surface source function [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn ! Radiances [W/m2-str]
      real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up ! Radiances [W/m2-str]
    end subroutine

    subroutine lw_source_2str(ncol, nlay, ngpt, top_at_1,   &
                              sfc_emis, sfc_src,      &
                              lay_source, lev_source, &
                              gamma1, gamma2, rdif, tdif, tau, source_dn, source_up, source_sfc) &
                              bind (C, name="lw_source_2str")
      use mo_rte_kind, only: wp, wl
      integer, value,                  intent(in) :: ncol, nlay, ngpt
      logical(wl), value,              intent(in) :: top_at_1
      real(wp), dimension(ncol      , ngpt), intent(in) :: sfc_emis, sfc_src
      real(wp), dimension(ncol, nlay, ngpt), intent(in) :: lay_source,    & ! Planck source at layer center
                                                     tau,           & ! Optical depth (tau)
                                                     gamma1, gamma2,& ! Coupling coefficients
                                                     rdif, tdif       ! Layer reflectance and transmittance
      real(wp), dimension(ncol, nlay+1, ngpt), target, &
                                       intent(in)  :: lev_source       ! Planck source at layer edges
      real(wp), dimension(ncol, nlay, ngpt), intent(out) :: source_dn, source_up
      real(wp), dimension(ncol      , ngpt), intent(out) :: source_sfc      ! Source function for upward radation at surface
    end subroutine

    subroutine lw_combine_sources(ncol, nlay, ngpt, top_at_1, &
                                  lev_src_inc, lev_src_dec, lev_source) bind(C, name="lw_combine_sources")
      use mo_rte_kind, only: wp, wl
      integer, value,                          intent(in ) :: ncol, nlay, ngpt
      logical(wl), value,                      intent(in ) :: top_at_1
      real(wp), dimension(ncol, nlay  , ngpt), intent(in ) :: lev_src_inc, lev_src_dec
      real(wp), dimension(ncol, nlay+1, ngpt), intent(out) :: lev_source
    end subroutine

    subroutine lw_two_stream(ncol, nlay, ngpt, tau, w0, g, &
                                  gamma1, gamma2, Rdif, Tdif) bind(C, name="lw_two_stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                      intent(in)  :: ncol, nlay, ngpt
      real(wp), dimension(ncol,nlay,ngpt), intent(in)  :: tau, w0, g
      real(wp), dimension(ncol,nlay,ngpt), intent(out) :: gamma1, gamma2, Rdif, Tdif
    end subroutine

    subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                                top_at_1, tau, mu0, flux_dir) bind (C, name="sw_solver_noscat")
      use mo_rte_kind, only: wp, wl
      integer, value,             intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,         intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
      real(wp), dimension(ncol            ), intent(in   ) :: mu0          ! cosine of solar zenith angle
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
    end subroutine

     subroutine lw_solver_2stream (ncol, nlay, ngpt, top_at_1, &
                                   tau, ssa, g,                &
                                   lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                                   flux_up, flux_dn) bind(C, name="lw_solver_2stream")
      use mo_rte_kind, only: wp, wl
      integer, value,                        intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl), value,                    intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau, &  ! Optical thickness,
                                                              ssa, &  ! single-scattering albedo,
                                                              g       ! asymmetry parameter []
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
      real(wp), dimension(ncol,nlay,ngpt), target, &
                                             intent(in   ) :: lev_source_inc, lev_source_dec
                                          ! Planck source at layer edge for radiation in increasing/decreasing ilay direction [W/m2]
                                          ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis         ! Surface emissivity      []
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src          ! Surface source function [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), &
                                             intent(  out) :: flux_up   ! Fluxes [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), &
                                             intent(inout) :: flux_dn  ! Fluxes [W/m2]
    end subroutine

  end interface

  interface apply_BC
    procedure apply_BC_0, apply_BC_factor, apply_BC_gpt
  end interface apply_BC

  public :: apply_BC, &
            lw_solver_noscat, lw_solver_noscat_GaussQuad, lw_solver_2stream, &
            sw_solver_noscat,                             sw_solver_2stream

  ! These routines don't really need to be visible but making them so is useful for testing.
  public :: lw_combine_sources, &
            lw_source_2str, sw_source_2str, &
            lw_two_stream, sw_two_stream, &
            adding

end module mo_rte_solver_kernels
