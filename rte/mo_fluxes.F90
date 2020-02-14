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
! Compute output quantities from RTE based on spectrally-resolved flux profiles
!    This module contains an class and a broadband implmentation that sums over all spectral points
!
! -------------------------------------------------------------------------------------------------
module mo_fluxes
  use mo_rte_kind,       only: wp
  use mo_rte_util_array, only: extents_are
  use mo_optical_props,  only: ty_optical_props
  use mo_fluxes_broadband_kernels, &
                         only: sum_broadband, net_broadband
  implicit none
  ! -----------------------------------------------------------------------------------------------
  !
  ! Base class
  !   reduce() function accepts spectral flux profiles, computes desired outputs
  !   are_desired() returns a logical - does it makes sense to invoke reduce()?
  !
  ! -----------------------------------------------------------------------------------------------
  type :: ty_fluxes
  end type ty_fluxes
  ! -----------------------------------------------------------------------------------------------
  !
  ! Class implementing broadband integration for the complete flux profile
  !   Data components are pointers so results can be written directly into memory
  !
  ! -----------------------------------------------------------------------------------------------
  type, extends(ty_fluxes) :: ty_fluxes_broadband
    real(wp), dimension(:,:), pointer :: flux_up => NULL(), flux_dn => NULL()
    real(wp), dimension(:,:), pointer :: flux_net => NULL()    ! Net (down - up)
    real(wp), dimension(:,:), pointer :: flux_dn_dir => NULL() ! Direct flux down
  end type ty_fluxes_broadband
  ! -----------------------------------------------------------------------------------------------
contains
  ! --------------------------------------------------------------------------------------
  !
  ! Broadband fluxes -- simply sum over the spectral dimension and report the whole profile
  !
  ! --------------------------------------------------------------------------------------
  function reduce(cls, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes_broadband),        intent(inout) :: cls
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt

    ! ------
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)
    error_msg = ""

    !
    ! Check array sizes
    !  Input arrays
    !
    if(.not. extents_are(gpt_flux_dn, ncol, nlev, ngpt)) then
      error_msg = "reduce: gpt_flux_dn array incorrectly sized"
      return
    end if
    if(present(gpt_flux_dn_dir)) then
      if(.not. extents_are(gpt_flux_dn_dir, ncol, nlev, ngpt)) then
        error_msg = "reduce: gpt_flux_dn_dir array incorrectly sized"
        return
      end if
    end if
    !
    ! Output arrays
    !
    if(associated(cls%flux_up)) then
      if(.not. extents_are(cls%flux_up, ncol, nlev)) then
        error_msg = 'reduce: flux_up array incorrectly sized'
        return
      end if
    end if
    if(associated(cls%flux_dn)) then
      if(.not. extents_are(cls%flux_dn, ncol, nlev)) then
        error_msg = 'reduce: flux_dn array incorrectly sized'
        return
      end if
    end if
    if(associated(cls%flux_net)) then
      if(.not. extents_are(cls%flux_net, ncol, nlev)) then
        error_msg = 'reduce: flux_net array incorrectly sized'
        return
      end if
    end if
    if(associated(cls%flux_dn_dir)) then
      if(.not. extents_are(cls%flux_dn_dir, ncol, nlev)) then
        error_msg = 'reduce: flux_dn_dir array incorrectly sized'
        return
      end if
    end if
    !
    ! Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if(associated(cls%flux_dn_dir) .and. .not. present(gpt_flux_dn_dir)) then
      error_msg = "reduce: requesting direct downward flux but cls hasn't been supplied"
      return
    end if

    !
    ! Broadband fluxes - call the kernels
    !
    if(associated(cls%flux_up    )) &
      call sum_broadband(ncol, nlev, ngpt, gpt_flux_up,     cls%flux_up)
    if(associated(cls%flux_dn    )) &
      call sum_broadband(ncol, nlev, ngpt, gpt_flux_dn,     cls%flux_dn)
    if(associated(cls%flux_dn_dir)) &
      call sum_broadband(ncol, nlev, ngpt, gpt_flux_dn_dir, cls%flux_dn_dir)

    if(associated(cls%flux_net)) then
      !
      !  Reuse down and up results if possible
      !
      if(associated(cls%flux_dn) .and. associated(cls%flux_up)) then
        call net_broadband(ncol, nlev,      cls%flux_dn, cls%flux_up, cls%flux_net)
      else
        call net_broadband(ncol, nlev, ngpt, gpt_flux_dn,  gpt_flux_up, cls%flux_net)
      end if
    end if
  end function reduce
  ! --------------------------------------------------------------------------------------
  !
  ! Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  !   be allocated for output
  !
  ! --------------------------------------------------------------------------------------
  function are_desired(cls)
    class(ty_fluxes_broadband), intent(in   ) :: cls
    logical                                   :: are_desired

    are_desired = any( [associated(cls%flux_up),     &
                        associated(cls%flux_dn),     &
                        associated(cls%flux_dn_dir), &
                        associated(cls%flux_net)] )
  end function are_desired
  ! --------------------------------------------------------------------------------------
end module mo_fluxes
