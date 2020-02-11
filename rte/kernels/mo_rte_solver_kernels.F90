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

  real(wp), parameter :: pi = acos(-1._wp)
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave two-stream calculation:
  !   combine RRTMGP-specific sources at levels
  !   compute layer reflectance, transmittance
  !   compute total source function at levels using linear-in-tau
  !   transport
  !
  ! -------------------------------------------------------------------------------------------------
   subroutine lw_solver_2stream (ncol, nlay, ngpt, top_at_1, &
                                 tau, ssa, g,                &
                                 lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                                 flux_up, flux_dn) bind(C, name="lw_solver_2stream")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
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
                                                                              ! Top level (= merge(1, nlay+1, top_at_1)
                                                                              ! must contain incident flux boundary condition
    ! ----------------------------------------------------------------------
    integer :: icol, igpt
    real(wp), dimension(ncol,nlay  ,ngpt) :: Rdif, Tdif, gamma1, gamma2
    real(wp), dimension(ncol       ,ngpt) :: sfc_albedo
    real(wp), dimension(ncol,nlay+1,ngpt) :: lev_source
    real(wp), dimension(ncol,nlay  ,ngpt) :: source_dn, source_up
    real(wp), dimension(ncol       ,ngpt) :: source_sfc
    ! ------------------------------------
    ! ------------------------------------
    !$acc enter data copyin(tau, ssa, g, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_dn)
    !$acc enter data create(flux_up, Rdif, Tdif, gamma1, gamma2, sfc_albedo, lev_source, source_dn, source_up, source_sfc)
    !
    ! RRTMGP provides source functions at each level using the spectral mapping
    !   of each adjacent layer. Combine these for two-stream calculations
    !
    call lw_combine_sources(ncol, nlay, ngpt, top_at_1, &
                            lev_source_inc, lev_source_dec, &
                            lev_source)
    !
    ! Cell properties: reflection, transmission for diffuse radiation
    !   Coupling coefficients needed for source function
    !
    call lw_two_stream(ncol, nlay, ngpt, &
                       tau , ssa, g,     &
                       gamma1, gamma2, Rdif, Tdif)

    !
    ! Source function for diffuse radiation
    !
    call lw_source_2str(ncol, nlay, ngpt, top_at_1, &
                        sfc_emis, sfc_src, &
                        lay_source, lev_source, &
                        gamma1, gamma2, Rdif, Tdif, tau, &
                        source_dn, source_up, source_sfc)

    !$acc  parallel loop collapse(2)
    do igpt = 1, ngpt
      do icol = 1, ncol
        sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt)
      end do
    end do
    !
    ! Transport
    !
    call adding(ncol, nlay, ngpt, top_at_1,        &
                sfc_albedo,                        &
                Rdif, Tdif,                        &
                source_dn, source_up, source_sfc,  &
                flux_up, flux_dn)
    !$acc exit data delete(tau, ssa, g, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src)
    !$acc exit data delete(Rdif, Tdif, gamma1, gamma2, sfc_albedo, lev_source, source_dn, source_up, source_sfc)
    !$acc exit data copyout(flux_up, flux_dn)
    write(*,*) 'GOT HERE: ', __FILE__, ": ", __LINE__
    write(*,*) "WARNING: THIS ISN'T TESTED"
  end subroutine lw_solver_2stream
  ! -------------------------------------------------------------------------------------------------
  !
  !   Top-level shortwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  !   Extinction-only i.e. solar direct beam
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                              top_at_1, tau, mu0, flux_dir) bind (C, name="sw_solver_noscat")
    integer,                    intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent(in   ) :: top_at_1
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol            ), intent(in   ) :: mu0          ! cosine of solar zenith angle
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
                                                                          ! Top level must contain incident flux boundary condition
    integer :: icol, ilev, igpt
    real(wp) :: mu0_inv(ncol)
    ! ------------------------------------
    ! ------------------------------------
    !$acc enter data copyin(tau, mu0) create(mu0_inv, flux_dir)
    !$acc parallel loop
    do icol = 1, ncol
      mu0_inv(icol) = 1._wp/mu0(icol)
    enddo
    ! Indexing into arrays for upward and downward propagation depends on the vertical
    !   orientation of the arrays (whether the domain top is at the first or last index)
    ! We write the loops out explicitly so compilers will have no trouble optimizing them.

    ! Downward propagation
    if(top_at_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      !$acc parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol
          do ilev = 2, nlay+1
            flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev-1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol))
          end do
        end do
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      !$acc parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol
          do ilev = nlay, 1, -1
            flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev+1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol))
          end do
        end do
      end do
    end if
    !$acc exit data delete(tau, mu0, mu0_inv) copyout(flux_dir)
    write(*,*) 'GOT HERE: ', __FILE__, ": ", __LINE__
    write(*,*) "WARNING: THIS ISN'T TESTED"
  end subroutine sw_solver_noscat

  ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
  !    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
  !
  ! Equations are developed in Meador and Weaver, 1980,
  !    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
  !
  subroutine lw_two_stream(ncol, nlay, ngpt, tau, w0, g, &
                                gamma1, gamma2, Rdif, Tdif) bind(C, name="lw_two_stream")
    integer,                             intent(in)  :: ncol, nlay, ngpt
    real(wp), dimension(ncol,nlay,ngpt), intent(in)  :: tau, w0, g
    real(wp), dimension(ncol,nlay,ngpt), intent(out) :: gamma1, gamma2, Rdif, Tdif

    ! -----------------------
    integer  :: icol, ilay, igpt

    ! Variables used in Meador and Weaver
    real(wp) :: k

    ! Ancillary variables
    real(wp) :: RT_term
    real(wp) :: exp_minusktau, exp_minus2ktau

    real(wp), parameter :: LW_diff_sec = 1.66  ! 1./cos(diffusivity angle)
    ! ---------------------------------
    ! ---------------------------------
    !$acc enter data copyin(tau, w0, g)
    !$acc enter data create(gamma1, gamma2, Rdif, Tdif)

    !$acc  parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          !
          ! Coefficients differ from SW implementation because the phase function is more isotropic
          !   Here we follow Fu et al. 1997, doi:10.1175/1520-0469(1997)054<2799:MSPITI>2.0.CO;2
          !   and use a diffusivity sec of 1.66
          !
          gamma1(icol,ilay,igpt)= LW_diff_sec * (1._wp - 0.5_wp * w0(icol,ilay,igpt) * (1._wp + g(icol,ilay,igpt))) ! Fu et al. Eq 2.9
          gamma2(icol,ilay,igpt)= LW_diff_sec *          0.5_wp * w0(icol,ilay,igpt) * (1._wp - g(icol,ilay,igpt))  ! Fu et al. Eq 2.10

          ! Written to encourage vectorization of exponential, square root
          ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
          !   k = 0 for isotropic, conservative scattering; this lower limit on k
          !   gives relative error with respect to conservative solution
          !   of < 0.1% in Rdif down to tau = 10^-9
          k = sqrt(max((gamma1(icol,ilay,igpt) - gamma2(icol,ilay,igpt)) * &
                       (gamma1(icol,ilay,igpt) + gamma2(icol,ilay,igpt)),  &
                       1.e-12_wp))
          exp_minusktau = exp(-tau(icol,ilay,igpt)*k)

          !
          ! Diffuse reflection and transmission
          !
          exp_minus2ktau = exp_minusktau * exp_minusktau

          ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
          RT_term = 1._wp / (k * (1._wp + exp_minus2ktau)  + &
                    gamma1(icol,ilay,igpt) * (1._wp - exp_minus2ktau) )

          ! Equation 25
          Rdif(icol,ilay,igpt) = RT_term * gamma2(icol,ilay,igpt) * (1._wp - exp_minus2ktau)

          ! Equation 26
          Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau
        end do
      end do
    end do
    !$acc exit data delete (tau, w0, g)
    !$acc exit data copyout(gamma1, gamma2, Rdif, Tdif)
    write(*,*) 'GOT HERE: ', __FILE__, ": ", __LINE__
    write(*,*) "WARNING: THIS ISN'T TESTED"
  end subroutine lw_two_stream
  ! -------------------------------------------------------------------------------------------------
  !
  ! Source function combination
  ! RRTMGP provides two source functions at each level
  !   using the spectral mapping from each of the adjascent layers.
  !   Need to combine these for use in two-stream calculation.
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_combine_sources(ncol, nlay, ngpt, top_at_1, &
                                lev_src_inc, lev_src_dec, lev_source) bind(C, name="lw_combine_sources")
    integer,                                 intent(in ) :: ncol, nlay, ngpt
    logical(wl),                             intent(in ) :: top_at_1
    real(wp), dimension(ncol, nlay  , ngpt), intent(in ) :: lev_src_inc, lev_src_dec
    real(wp), dimension(ncol, nlay+1, ngpt), intent(out) :: lev_source

    integer :: icol, ilay, igpt
    ! ---------------------------------------------------------------
    ! ---------------------------------
    !$acc enter data copyin(lev_src_inc, lev_src_dec)
    !$acc enter data create(lev_source)

    !$acc  parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay+1
        do icol = 1,ncol
          if(ilay == 1) then
            lev_source(icol, ilay, igpt) =      lev_src_dec(icol, ilay,   igpt)
          else if (ilay == nlay+1) then
            lev_source(icol, ilay, igpt) =      lev_src_inc(icol, ilay-1, igpt)
          else
            lev_source(icol, ilay, igpt) = sqrt(lev_src_dec(icol, ilay, igpt) * &
                                                lev_src_inc(icol, ilay-1, igpt))
          end if
        end do
      end do
    end do
    !$acc exit data delete (lev_src_inc, lev_src_dec)
    !$acc exit data copyout(lev_source)
    write(*,*) 'GOT HERE: ', __FILE__, ": ", __LINE__
    write(*,*) "WARNING: THIS ISN'T TESTED"
  end subroutine lw_combine_sources
  ! ---------------------------------------------------------------
  !
  ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
  !   This version straight from ECRAD
  !   Source is provided as W/m2-str; factor of pi converts to flux units
  !
  ! ---------------------------------------------------------------
  subroutine lw_source_2str(ncol, nlay, ngpt, top_at_1,   &
                            sfc_emis, sfc_src,      &
                            lay_source, lev_source, &
                            gamma1, gamma2, rdif, tdif, tau, source_dn, source_up, source_sfc) &
                            bind (C, name="lw_source_2str")
    integer,                         intent(in) :: ncol, nlay, ngpt
    logical(wl),                     intent(in) :: top_at_1
    real(wp), dimension(ncol      , ngpt), intent(in) :: sfc_emis, sfc_src
    real(wp), dimension(ncol, nlay, ngpt), intent(in) :: lay_source,    & ! Planck source at layer center
                                                   tau,           & ! Optical depth (tau)
                                                   gamma1, gamma2,& ! Coupling coefficients
                                                   rdif, tdif       ! Layer reflectance and transmittance
    real(wp), dimension(ncol, nlay+1, ngpt), target, &
                                     intent(in)  :: lev_source       ! Planck source at layer edges
    real(wp), dimension(ncol, nlay, ngpt), intent(out) :: source_dn, source_up
    real(wp), dimension(ncol      , ngpt), intent(out) :: source_sfc      ! Source function for upward radation at surface

    integer             :: icol, ilay, igpt
    real(wp)            :: Z, Zup_top, Zup_bottom, Zdn_top, Zdn_bottom
    real(wp)            :: lev_source_bot, lev_source_top
    ! ---------------------------------------------------------------
    ! ---------------------------------
    !$acc enter data copyin(sfc_emis, sfc_src, lay_source, tau, gamma1, gamma2, rdif, tdif, lev_source)
    !$acc enter data create(source_dn, source_up, source_sfc)

    !$acc parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          if (tau(icol,ilay,ngpt) > 1.0e-8_wp) then
            if(top_at_1) then
              lev_source_top = lev_source(icol,ilay  ,ngpt)
              lev_source_bot = lev_source(icol,ilay+1,ngpt)
            else
              lev_source_top = lev_source(icol,ilay+1,ngpt)
              lev_source_bot = lev_source(icol,ilay  ,ngpt)
            end if
            !
            ! Toon et al. (JGR 1989) Eqs 26-27
            !
            Z = (lev_source_bot-lev_source_top) / (tau(icol,ilay,igpt)*(gamma1(icol,ilay,igpt)+gamma2(icol,ilay,igpt)))
            Zup_top        =  Z + lev_source_top
            Zup_bottom     =  Z + lev_source_bot
            Zdn_top        = -Z + lev_source_top
            Zdn_bottom     = -Z + lev_source_bot
            source_up(icol,ilay,igpt) = pi * (Zup_top    - rdif(icol,ilay,igpt) * Zdn_top    - tdif(icol,ilay,igpt) * Zup_bottom)
            source_dn(icol,ilay,igpt) = pi * (Zdn_bottom - rdif(icol,ilay,igpt) * Zup_bottom - tdif(icol,ilay,igpt) * Zdn_top)
          else
            source_up(icol,ilay,igpt) = 0._wp
            source_dn(icol,ilay,igpt) = 0._wp
          end if
          if(ilay == 1) source_sfc(icol,igpt) = pi * sfc_emis(icol,igpt) * sfc_src(icol,igpt)
        end do
      end do
    end do
    !$acc exit data delete(sfc_emis, sfc_src, lay_source, tau, gamma1, gamma2, rdif, tdif, lev_source)
    !$acc exit data copyout(source_dn, source_up, source_sfc)
    write(*,*) 'GOT HERE: ', __FILE__, ": ", __LINE__
    write(*,*) "WARNING: THIS ISN'T TESTED"
  end subroutine lw_source_2str

end module mo_rte_solver_kernels
