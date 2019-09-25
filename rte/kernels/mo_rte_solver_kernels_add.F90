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
module mo_rte_solver_kernels_add
  use,  intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  use mo_rte_solver_kernels, only : apply_BC
  implicit none
  private

  public :: lw_solver_Tang_GaussQuad, lw_solver_Tang
  public :: lw_solver_Tip_GaussQuad,  lw_solver_Tip

  real(wp), parameter :: pi = acos(-1._wp)
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! Top-level longwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW fluxes, no scattering, mu (cosine of integration angle) specified by column
  !   Does radiation calculation at user-supplied angles; converts radiances to flux
  !   using user-supplied weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_Tang(ncol, nlay, ngpt, top_at_1, D, weight,                             &
                              tau, ssa, lay_source, lev_source_inc, &
                              lev_source_dec, sfc_emis, sfc_src, &
                              planckLev, planckSfc,plancFrac,&
                              radn_up, radn_dn) &
                              bind(C, name="lw_solver_Tang")
  use mo_rte_solver_kernels, only : lw_source_noscat
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: D            ! secant of propagation angle  []
    real(wp),                              intent(in   ) :: weight       ! quadrature weight
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: ssa          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
    ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
    ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
    real(wp), dimension(ncol,nlay,  ngpt), target, &
                                           intent(in   ) :: lev_source_inc, lev_source_dec
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: planckLev        ! Level Plank source      [w/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: planckSfc        ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: plancFrac        ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn      ! Top level must contain incident flux boundary condition

    ! Local variables, no g-point dependency
    real(wp), dimension(ncol,nlay) :: tau_loc, &  ! path length (tau/mu)
                                        trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay) :: source_dn, source_up
    real(wp), dimension(ncol     ) :: source_sfc, sfc_albedo

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    real(wp), parameter :: pi = acos(-1._wp)
    integer             :: ilev, igpt, top_level
    ! ------------------------------------
    real(wp)            :: fact
    real(wp), parameter :: tau_thresh = sqrt(epsilon(tau))
    integer             :: icol
    real(wp), dimension(ncol     ) :: sfcSource
    ! ------------------------------------

    ! Which way is up?
    ! Level Planck sources for upward and downward radiation
    ! When top_at_1, lev_source_up => lev_source_dec
    !                lev_source_dn => lev_source_inc, and vice-versa
    if(top_at_1) then
      top_level = 1
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      top_level = nlay+1
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    do igpt = 1, ngpt
      !
      ! Transport is for intensity
      !   convert flux at top of domain to intensity assuming azimuthal isotropy
      !
      radn_dn(:,top_level,igpt) = radn_dn(:,top_level,igpt)/(2._wp * pi * weight)

      !
      ! Optical path and transmission, used in source function and transport calculations
      !
      do ilev = 1, nlay
        tau_loc(:,ilev) = tau(:,ilev,igpt)*D(:,igpt)
        trans  (:,ilev) = exp(-tau_loc(:,ilev))
      end do

      !
      ! Source function for diffuse radiation
      !
      ! here lev_source_up and lev_source_dn are level source multiplied by 
      ! layer planck fraction
      !
      call lw_source_noscat(ncol, nlay, lay_source(:,:,igpt), &
                            lev_source_up(:,:,igpt), lev_source_dn(:,:,igpt), &
                            tau_loc, trans, source_dn, source_up)
      
      !
      ! Surface albedo, surface source function
      !
      sfc_albedo(:) = 1._wp - sfc_emis(:,igpt)
      source_sfc(:) = sfc_emis(:,igpt) * sfc_src(:,igpt)
      sfcSource     = planckSfc(:,igpt) * sfc_emis(:,igpt)
      !
      ! Transport
      !
      call lw_transport_Tang(ncol, nlay, top_at_1,  &
                             tau_loc, ssa(:,:,igpt), trans, &
                             sfc_albedo, source_dn, source_up, source_sfc, &
                             radn_up(:,:,igpt), radn_dn(:,:,igpt), &
                             planckLev(:,:,igpt),&
                             plancFrac(:,:,igpt))

     
      ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      radn_dn(:,:,igpt) = 2._wp * pi * weight * radn_dn(:,:,igpt)
      radn_up(:,:,igpt) = 2._wp * pi * weight * radn_up(:,:,igpt)
    end do  ! g point loop

  end subroutine lw_solver_Tang
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !   Routine sums over single-angle solutions for each sets of angles/weights
  !
  ! ---------------------------------------------------------------
 
  subroutine lw_solver_Tang_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                   tau, ssa, lay_source, lev_source_inc, lev_source_dec, &
                                   sfc_emis, sfc_src,planckLev, planckSfc,&
                                   plancFrac,&
                                  flux_up, flux_dn) &
                                   bind(C, name="lw_solver_Tang_GaussQuad")
    use mo_rte_solver_kernels, only : lw_solver_noscat
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    integer,                               intent(in   ) :: nmus         ! number of quadrature angles
    real(wp), dimension(nmus),             intent(in   ) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: ssa          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: planckLev        ! Level Plank source      [w/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: planckSfc        ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: plancFrac        ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn      ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol

    integer :: imu, top_level
    ! ------------------------------------
    !
    ! For the first angle output arrays store total flux
    !
    Ds_ncol(:,:) = Ds(1)
    call lw_solver_Tang(ncol, nlay, ngpt, &
                          top_at_1, Ds_ncol, weights(1), tau, ssa, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          planckLev, planckSfc, plancFrac,flux_up, flux_dn)

    !
    ! For more than one angle use local arrays
    !
    top_level = MERGE(1, nlay+1, top_at_1)
    call apply_BC(ncol, nlay, ngpt, top_at_1, flux_dn(:,top_level,:), radn_dn)

    do imu = 2, nmus
      Ds_ncol(:,:) = Ds(imu)
      call lw_solver_Tang(ncol, nlay, ngpt, &
                            top_at_1, Ds_ncol, weights(imu), tau, ssa, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            planckLev, planckSfc, plancFrac, radn_up, radn_dn)

      flux_up(:,:,:) = flux_up(:,:,:) + radn_up(:,:,:)
      flux_dn(:,:,:) = flux_dn(:,:,:) + radn_dn(:,:,:)
    end do
  end subroutine lw_solver_Tang_GaussQuad

    ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave no-scattering transport
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_transport_Tang(ncol, nlay, top_at_1, &
                                 tau, ssa, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                 radn_up, radn_dn,&
                                 planckLev, planckFrc) bind(C, name="lw_transport_Tang")
    integer,                          intent(in   ) :: ncol, nlay ! Number of columns, layers, g-points
    logical(wl),                      intent(in   ) :: top_at_1   !
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                       trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: ssa        ! single scattering
    real(wp), dimension(ncol       ), intent(in   ) :: sfc_albedo ! Surface albedo
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: source_dn, &
                                                       source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ), intent(in   ) :: source_sfc ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1), intent(in   ) :: planckLev  ! Level Plank source      [w/m2]
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: planckFrc  ! Layer Plank fraction    [1]
    real(wp), dimension(ncol,nlay+1), intent(inout) :: radn_up    ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1), intent(inout) :: radn_dn    !Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev, icol
    real(wp), dimension(ncol,nlay+1) :: i_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1) :: i_dn      ! Top level must contain incident flux boundary condition

    real(wp)   :: xxx, yyy
    ! ---------------------------------------------------
    if(top_at_1) then
      !
      ! Top of domain is index 1
      !
      ! Downward propagation
      i_dn(:,1)=radn_dn(:,1)/planckFrc(:,1)

      do ilev = 1, nlay
        radn_dn(:,ilev+1) = trans(:,ilev)*radn_dn(:,ilev) + source_dn(:,ilev)
        i_dn(:,ilev+1)    = radn_dn(:,ilev+1)/planckFrc(:,ilev)
      end do

      ! Surface reflection and emission
      radn_up(:,nlay+1) = radn_dn(:,nlay+1)*sfc_albedo(:) + source_sfc(:)
      i_up   (:,nlay+1) = radn_up(:,nlay+1)/planckFrc(:,nlay)

      ! Upward propagation
      do ilev = nlay, 1, -1
        radn_up(:,ilev) = trans(:,ilev  )*radn_up(:,ilev+1) + source_up(:,ilev)
        i_up(:,ilev)    = radn_up(:,ilev)/planckFrc(:,ilev)
        do icol=1,ncol
           if ( ssa(icol,ilev) > 1e-6 )  then
              radn_up(icol,ilev) = radn_up(icol,ilev) + &
                      0.5_wp*ssa(icol,ilev)*planckFrc(icol,ilev) *&
                      ( i_dn(icol,ilev)   - planckLev(icol,ilev) -  &
                       (i_dn(icol,ilev+1) - planckLev(icol,ilev+1))*trans(icol,ilev) )
          endif  
        enddo  
      end do
      ! 2nd Downward propagation
      do ilev = 1, nlay
        radn_dn(:,ilev+1) = trans(:,ilev)*radn_dn(:,ilev) + source_dn(:,ilev)
        do icol=1,ncol
          if ( ssa(icol,ilev) > 1e-6 )  then
            radn_dn(icol,ilev) = radn_dn(icol,ilev) + &
                   0.5_wp*ssa(icol,ilev)*planckFrc(icol,ilev) * &
                   ( i_up(icol,ilev)   - planckLev(icol,ilev) -  &
                    (i_up(icol,ilev+1) - planckLev(icol,ilev+1))*trans(icol,ilev) )
          endif  
        end do
      end do
    else
      !
      ! Top of domain is index nlay+1
      !
      ! Downward propagation
      i_dn(:,nlay+1)=radn_dn(:,nlay+1)/planckFrc(:,nlay)
      do ilev = nlay, 1, -1
        radn_dn(:,ilev) = trans(:,ilev  )*radn_dn(:,ilev+1) + source_dn(:,ilev)
        i_dn(:,ilev)    = radn_dn(:,ilev)/planckFrc(:,ilev)
      end do

      ! Surface reflection and emission
      radn_up(:,     1) = radn_dn(:,     1)*sfc_albedo(:) + source_sfc(:)
      i_up   (:,     1) = radn_dn(:,     1)/planckFrc(:,1)

      ! Upward propagation
      do ilev = 1, nlay
        radn_up(:,ilev+1) = trans(:,ilev) * radn_up(:,ilev) +  source_up(:,ilev)
        i_up   (:,ilev+1) = radn_up(:,ilev+1)/planckFrc(:,ilev)
        do icol=1,ncol
           if ( ssa(icol,ilev) > 1e-6 )  then
              yyy = ( i_dn(icol,ilev+1)   - planckLev(icol,ilev+1) -  &
                      (i_dn(icol,ilev) - planckLev(icol,ilev))*trans(icol,ilev) )

              xxx = 0.5_wp*ssa(icol,ilev)*planckFrc(icol,ilev) * yyy
              radn_up(icol,ilev+1) = radn_up(icol,ilev+1) + xxx
           endif  
        enddo  
      end do

      ! 2st Downward propagation
      do ilev = nlay, 1, -1
        radn_dn(:,ilev) = trans(:,ilev  )*radn_dn(:,ilev+1) + source_dn(:,ilev)
        do icol=1,ncol
          if ( ssa(icol,ilev) > 1e-6 )  then
            radn_dn(icol,ilev) = radn_dn(icol,ilev) + &
                   0.5_wp*ssa(icol,ilev)*planckFrc(icol,ilev) * &
                   ( i_up(icol,ilev)   - planckLev(icol,ilev) -  &
                    (i_up(icol,ilev+1) - planckLev(icol,ilev+1))*trans(icol,ilev) )
          endif  
        enddo  
      end do
    end if
  end subroutine lw_transport_Tang

 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  TANG + IP
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! -------------------------------------------------------------------------------------------------
  !
  ! Top-level longwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW fluxes, no scattering, mu (cosine of integration angle) specified by column
  !   Does radiation calculation at user-supplied angles; converts radiances to flux
  !   using user-supplied weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_Tip(ncol, nlay, ngpt, top_at_1, D, weight,                             &
                              tau, ssa, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                              radn_up, radn_dn) bind(C, name="lw_solver_Tip")
  use mo_rte_solver_kernels, only : lw_source_noscat, lw_transport_noscat

    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: D            ! secant of propagation angle  []
    real(wp),                              intent(in   ) :: weight       ! quadrature weight
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: ssa          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
    ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
    ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
    real(wp), dimension(ncol,nlay,  ngpt), target, &
                                           intent(in   ) :: lev_source_inc, lev_source_dec
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn      ! Top level must contain incident flux boundary condition

    ! Local variables, no g-point dependency
    real(wp), dimension(ncol,nlay) :: tau_loc, &  ! path length (tau/mu)
                                        trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay) :: source_dn, source_up
    real(wp), dimension(ncol     ) :: source_sfc, sfc_albedo

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    real(wp), parameter :: pi = acos(-1._wp)
    integer             :: ilev, igpt, top_level
    ! ------------------------------------
    real(wp)            :: fact
    real(wp), parameter :: tau_thresh = sqrt(epsilon(tau))
    integer             :: icol
    real(wp), dimension(ncol     ) :: sfcSource
    ! ------------------------------------

    ! Which way is up?
    ! Level Planck sources for upward and downward radiation
    ! When top_at_1, lev_source_up => lev_source_dec
    !                lev_source_dn => lev_source_inc, and vice-versa
    if(top_at_1) then
      top_level = 1
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      top_level = nlay+1
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    do igpt = 1, ngpt
      !
      ! Transport is for intensity
      !   convert flux at top of domain to intensity assuming azimuthal isotropy
      !
      radn_dn(:,top_level,igpt) = radn_dn(:,top_level,igpt)/(2._wp * pi * weight)

      !
      ! Optical path and transmission, used in source function and transport calculations
      !
      !
      ! Optical path and transmission, used in source function and transport calculations
      !
      do ilev = 1, nlay
        tau_loc(:,ilev) = tau(:,ilev,igpt)*D(:,igpt)
        trans  (:,ilev) = exp(-tau_loc(:,ilev))
      end do

      ! forall(icol = 1:ncol,ilev = 1:nlay) tau_loc(icol,ilev) = tau(icol,ilev,igpt)*D(icol,igpt)
      ! forall(icol = 1:ncol,ilev = 1:nlay) trans  (icol,ilev) = exp(-tau_loc(icol,ilev))
      !
      ! Source function for diffuse radiation
      !
      call lw_source_noscat(ncol, nlay, &
                            lay_source(:,:,igpt), lev_source_up(:,:,igpt), lev_source_dn(:,:,igpt), &
                            tau_loc, trans, source_dn, source_up)

      !
      ! Surface albedo, surface source function
      !
      sfc_albedo(:) = 1._wp - sfc_emis(:,igpt)
      source_sfc(:) = sfc_emis(:,igpt) * sfc_src(:,igpt)
      !
      ! Transport
      !
      call lw_transport_Tip(ncol, nlay, top_at_1,  &
                               tau_loc, ssa(:,:,igpt), trans, &
                               sfc_albedo, source_dn, source_up, source_sfc, &
                               radn_up(:,:,igpt), radn_dn(:,:,igpt))

      ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      radn_dn(:,:,igpt) = 2._wp * pi * weight * radn_dn(:,:,igpt)
      radn_up(:,:,igpt) = 2._wp * pi * weight * radn_up(:,:,igpt)
    end do  ! g point loop
  end subroutine lw_solver_Tip
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !   Routine sums over single-angle solutions for each sets of angles/weights
  !
  ! ---------------------------------------------------------------
 
  subroutine lw_solver_Tip_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                   tau, ssa, lay_source, lev_source_inc, lev_source_dec, &
                                   sfc_emis, sfc_src,&
                                  flux_up, flux_dn) &
                                   bind(C, name="lw_solver_Tip_GaussQuad")
    use mo_rte_solver_kernels, only : lw_solver_noscat
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    integer,                               intent(in   ) :: nmus         ! number of quadrature angles
    real(wp), dimension(nmus),             intent(in   ) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: ssa          ! single scattering albedo []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay+1,ngpt), intent(in   ) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis     ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src      ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up      ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn      ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol

    integer :: imu, top_level
    ! ------------------------------------
    !
    ! For the first angle output arrays store total flux
    !
    Ds_ncol(:,:) = Ds(1)
    call lw_solver_Tip(ncol, nlay, ngpt, &
                          top_at_1, Ds_ncol, weights(1), tau, ssa, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          flux_up, flux_dn)

    !
    ! For more than one angle use local arrays
    !
    top_level = MERGE(1, nlay+1, top_at_1)
    ! call apply_BC(ncol, nlay, ngpt, top_at_1, flux_dn(:,top_level,:), radn_dn)
    radn_dn(1:ncol, top_level, 1:ngpt)  = flux_dn(1:ncol, top_level, 1:ngpt)

    do imu = 2, nmus
      Ds_ncol(:,:) = Ds(imu)
      call lw_solver_Tip(ncol, nlay, ngpt, &
                            top_at_1, Ds_ncol, weights(imu), tau, ssa, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            radn_up, radn_dn)

      flux_up(:,:,:) = flux_up(:,:,:) + radn_up(:,:,:)
      flux_dn(:,:,:) = flux_dn(:,:,:) + radn_dn(:,:,:)
    end do
  end subroutine lw_solver_Tip_GaussQuad

    ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave no-scattering transport
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_transport_Tip(ncol, nlay, top_at_1, &
                                 tau, ssa, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                 radn_up, radn_dn) bind(C, name="lw_transport_Tip")
    integer,                          intent(in   ) :: ncol, nlay ! Number of columns, layers, g-points
    logical(wl),                      intent(in   ) :: top_at_1   !
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                       trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: ssa        ! single scattering
    real(wp), dimension(ncol       ), intent(in   ) :: sfc_albedo ! Surface albedo
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: source_dn, &
                                                       source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ), intent(in   ) :: source_sfc ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1), intent(inout) :: radn_up    ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1), intent(inout) :: radn_dn    !Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev, icol
    ! ---------------------------------------------------
    real(wp), dimension(ncol,nlay+1)                :: i_up    ! Radiances [W/m2-str]
    real(wp) :: xx, yy, xx1, xx2, xx3
    if(top_at_1) then
      !
      ! Top of domain is index 1
      !
      ! Downward propagation
      do ilev = 1, nlay
        radn_dn(:,ilev+1) = trans(:,ilev)  *radn_dn(:,ilev)  + source_dn(:,ilev)
      end do

      ! Surface reflection and emission
      radn_up(:,nlay+1) = radn_dn(:,nlay+1)*sfc_albedo(:) + source_sfc(:)

      ! Surface reflection and emission
      radn_up(:,     1) = radn_dn(:,     1)*sfc_albedo(:) + source_sfc(:)
      i_up   (:,     1) = radn_up(:,     1)

      ! 1st Upward propagation
      do ilev = nlay, 1, -1
        radn_up(:,ilev) = trans(:,ilev  )*radn_up(:,ilev+1) + source_up(:,ilev)
        i_up   (:,ilev) = radn_up(:,ilev)

        do icol=1,ncol
           if ( ssa(icol,ilev) > 1e-6 )  then
              xx = 0.5_wp*ssa(icol,ilev)*&
                     ( radn_dn(icol,ilev)*(1.-trans(icol,ilev )**2 ) - &
                       source_dn(icol,ilev) *trans(icol,ilev  ) - &
                       source_up(icol,ilev))
              radn_up(icol,ilev) = radn_up(icol,ilev) + xx
            endif  
          enddo  
      end do
      ! 2nd Downward propagation
      do ilev = 1, nlay
        radn_dn(:,ilev+1) = trans(:,ilev)*radn_dn(:,ilev) + source_dn(:,ilev)
        do icol=1,ncol
          if ( ssa(icol,ilev) > 1e-6 )  then
              xx = 0.5_wp*ssa(icol,ilev)*( &
                  i_up(icol,ilev)*(1.-trans(icol,ilev  )**2)  - &
                  source_up(icol,ilev)*trans(icol,ilev  ) - &
                  source_dn(icol,ilev) )
                radn_dn(:,ilev+1) = radn_dn(:,ilev+1) + xx
          endif  
        enddo  
      end do
    else
      !
      ! Top of domain is index nlay+1
      !
      ! Downward propagation
      !
      ! --------- N+1
      !                   layer N
      ! ----------N
      !
      !
      do ilev = nlay, 1, -1
        radn_dn(:,ilev) = trans(:,ilev  )*radn_dn(:,ilev+1) + source_dn(:,ilev)
      end do

      ! Surface reflection and emission
      radn_up(:,     1) = radn_dn(:,     1)*sfc_albedo(:) + source_sfc(:)
      i_up   (:,     1) = radn_up(:,     1)
      ! Upward propagation
      do ilev = 1, nlay
        radn_up(:,ilev+1) =  trans(:,ilev) * radn_up(:,ilev) +  source_up(:,ilev)
        i_up   (:,ilev+1) = radn_up(:,ilev+1)
        do icol=1,ncol
           if ( ssa(icol,ilev) > 1e-6 )  then
              xx = 0.5_wp*ssa(icol,ilev)*&
                     ( radn_dn(icol,ilev+1)*(1.-trans(icol,ilev )**2 ) - &
                       source_dn(icol,ilev) *trans(icol,ilev  ) - &
                       source_up(icol,ilev))
              radn_up(icol,ilev+1) = radn_up(icol,ilev+1) + xx
          endif  
        enddo  
      end do

      ! 2st Downward propagation
      do ilev = nlay, 1, -1
        radn_dn(:,ilev) = trans(:,ilev  )*radn_dn(:,ilev+1) + source_dn(:,ilev)
        do icol=1,ncol
           if ( ssa(icol,ilev) > 1e-6 )  then
                     xx = 0.5_wp*ssa(icol,ilev)*( &
                      i_up(icol,ilev)*(1.-trans(icol,ilev  )**2)  - &
                      source_up(icol,ilev)*trans(icol,ilev  ) - &
                      source_dn(icol,ilev) )
              radn_dn(icol,ilev) = radn_dn(icol,ilev) + xx
          endif  
        enddo  
      end do
    end if
  end subroutine lw_transport_Tip
end module mo_rte_solver_kernels_add
