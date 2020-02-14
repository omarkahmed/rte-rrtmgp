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
! Encapsulate source function arrays for longwave/lw/internal sources
!    and shortwave/sw/external source.
!
! -------------------------------------------------------------------------------------------------
module mo_source_functions
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props, get_ngpt, init, is_initialized, finalize_props => finalize
  implicit none
  ! -------------------------------------------------------------------------------------------------
  !
  ! Type for longwave sources: computed at layer center, at layer edges using
  !   spectral mapping in each direction separately, and at the surface
  !
  type, extends(ty_optical_props), public :: ty_source_func_lw
    real(wp), allocatable, dimension(:,:,:) :: lay_source,     & ! Planck source at layer average temperature
                                                                 ! [W/m2] (ncol, nlay, ngpt)
                                               lev_source_inc, &  ! Planck source at layer edge,
                                               lev_source_dec     ! [W/m2] (ncol, nlay, ngpt)
                                                                  ! in increasing/decreasing ilay direction
                                                                  ! Includes spectral weighting that accounts for state-dependent
                                                                  ! frequency to g-space mapping
    real(wp), allocatable, dimension(:,:  ) :: sfc_source
  end type ty_source_func_lw
  ! -------------------------------------------------------------------------------------------------
  !
  ! Type for shortave sources: top-of-domain spectrally-resolved flux
  !
  type, extends(ty_optical_props), public :: ty_source_func_sw
    real(wp), allocatable, dimension(:,:  ) :: toa_source
  end type ty_source_func_sw

  interface finalize
    module procedure :: finalize_lw, finalize_sw
  end interface finalize

  interface is_allocated
    module procedure :: is_allocated_lw, is_allocated_sw
  end interface

  interface get_subset
    module procedure :: get_subset_range_sw, get_subset_range_lw
  end interface

  interface get_ncol
    module procedure :: get_ncol_lw, get_ncol_sw
  end interface

  interface alloc
    module procedure :: alloc_lw, copy_and_alloc_lw, alloc_sw, copy_and_alloc_sw
  end interface alloc


  ! -------------------------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for initialization, validity checking, finalization
  !
  ! ------------------------------------------------------------------------------------------
  !
  ! Longwave
  !
  ! ------------------------------------------------------------------------------------------
  pure function is_allocated_lw(cls)
    class(ty_source_func_lw), intent(in) :: cls
    logical                              :: is_allocated_lw

    is_allocated_lw = is_initialized(cls) .and. &
                      allocated(cls%sfc_source)
  end function is_allocated_lw
  ! --------------------------------------------------------------
  function alloc_lw(cls, ncol, nlay) result(err_message)
    class(ty_source_func_lw),    intent(inout) :: cls
    integer,                     intent(in   ) :: ncol, nlay
    character(len = 128)                       :: err_message

    integer :: ngpt

    err_message = ""
    if(.not. is_initialized(cls)) &
      err_message = "source_func_lw%alloc: not initialized so can't allocate"
    if(any([ncol, nlay] <= 0)) &
      err_message = "source_func_lw%alloc: must provide positive extents for ncol, nlay"
    if (err_message /= "") return

    if(allocated(cls%sfc_source)) deallocate(cls%sfc_source)
    if(allocated(cls%lay_source)) deallocate(cls%lay_source)
    if(allocated(cls%lev_source_inc)) deallocate(cls%lev_source_inc)
    if(allocated(cls%lev_source_dec)) deallocate(cls%lev_source_dec)

    ngpt = get_ngpt(cls)
    allocate(cls%sfc_source    (ncol,     ngpt), cls%lay_source    (ncol,nlay,ngpt), &
             cls%lev_source_inc(ncol,nlay,ngpt), cls%lev_source_dec(ncol,nlay,ngpt))
  end function alloc_lw
  ! --------------------------------------------------------------
  function copy_and_alloc_lw(cls, ncol, nlay, spectral_desc) result(err_message)
    class(ty_source_func_lw),    intent(inout) :: cls
    integer,                     intent(in   ) :: ncol, nlay
    class(ty_optical_props ),    intent(in   ) :: spectral_desc
    character(len = 128)                       :: err_message

    err_message = ""
    if(.not. is_initialized(spectral_desc)) then
      err_message = "source_func_lw%alloc: spectral_desc not initialized"
      return
    end if
    call finalize(cls)
    err_message = init(cls,spectral_desc)
    if (err_message /= "") return
    err_message = alloc(cls,ncol,nlay)
  end function copy_and_alloc_lw
  ! ------------------------------------------------------------------------------------------
  !
  ! Shortwave
  !
  ! ------------------------------------------------------------------------------------------
  pure function is_allocated_sw(cls)
    class(ty_source_func_sw), intent(in) :: cls
    logical                              :: is_allocated_sw

    is_allocated_sw = is_initialized(cls%ty_optical_props) .and. &
                      allocated(cls%toa_source)
  end function is_allocated_sw
  ! --------------------------------------------------------------
  function alloc_sw(cls, ncol) result(err_message)
    class(ty_source_func_sw),    intent(inout) :: cls
    integer,                     intent(in   ) :: ncol
    character(len = 128)                       :: err_message

    err_message = ""
    if(.not. is_initialized(cls)) &
      err_message = "source_func_sw%alloc: not initialized so can't allocate"
    if(ncol <= 0) &
      err_message = "source_func_sw%alloc: must provide positive extents for ncol"
    if (err_message /= "") return

    if(allocated(cls%toa_source)) deallocate(cls%toa_source)

    allocate(cls%toa_source(ncol, get_ngpt(cls)))
  end function alloc_sw
  ! --------------------------------------------------------------
  function copy_and_alloc_sw(cls, ncol, spectral_desc) result(err_message)
    class(ty_source_func_sw),    intent(inout) :: cls
    integer,                     intent(in   ) :: ncol
    class(ty_optical_props ),    intent(in   ) :: spectral_desc
    character(len = 128)                       :: err_message

    err_message = ""
    if(.not. is_initialized(spectral_desc)) then
      err_message = "source_func_sw%alloc: spectral_desc not initialized"
      return
    end if
    err_message = init(cls,spectral_desc)
    if(err_message /= "") return
    err_message = alloc(cls,ncol)
  end function copy_and_alloc_sw
  ! ------------------------------------------------------------------------------------------
  !
  ! Finalization (memory deallocation)
  !
  ! ------------------------------------------------------------------------------------------
  subroutine finalize_lw(cls)
    class(ty_source_func_lw),    intent(inout) :: cls

    if(allocated(cls%lay_source    )) deallocate(cls%lay_source)
    if(allocated(cls%lev_source_inc)) deallocate(cls%lev_source_inc)
    if(allocated(cls%lev_source_dec)) deallocate(cls%lev_source_dec)
    if(allocated(cls%sfc_source    )) deallocate(cls%sfc_source)
    call finalize_props(cls%ty_optical_props)
  end subroutine finalize_lw
  ! --------------------------------------------------------------
  subroutine finalize_sw(cls)
    class(ty_source_func_sw),    intent(inout) :: cls

    if(allocated(cls%toa_source    )) deallocate(cls%toa_source)
    call finalize_props(cls%ty_optical_props)
  end subroutine finalize_sw
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for finding the problem size
  !
  ! ------------------------------------------------------------------------------------------
  pure function get_ncol_lw(cls)
    class(ty_source_func_lw), intent(in) :: cls
    integer :: get_ncol_lw

    if(is_allocated(cls)) then
      get_ncol_lw = size(cls%lay_source,1)
    else
      get_ncol_lw = 0
    end if
  end function get_ncol_lw
  ! --------------------------------------------------------------
  pure function get_nlay(cls)
    class(ty_source_func_lw), intent(in) :: cls
    integer :: get_nlay

    if(is_allocated(cls)) then
      get_nlay = size(cls%lay_source,2)
    else
      get_nlay = 0
    end if
  end function get_nlay
  ! --------------------------------------------------------------
  pure function get_ncol_sw(cls)
    class(ty_source_func_sw), intent(in) :: cls
    integer :: get_ncol_sw

    if(is_allocated(cls)) then
      get_ncol_sw = size(cls%toa_source,1)
    else
      get_ncol_sw = 0
    end if
  end function get_ncol_sw
  ! ------------------------------------------------------------------------------------------
  !
  !  Routines for subsetting
  !
  ! ------------------------------------------------------------------------------------------
  function get_subset_range_lw(full, start, n, subset) result(err_message)
    class(ty_source_func_lw), intent(inout) :: full
    integer,                  intent(in   ) :: start, n
    class(ty_source_func_lw), intent(inout) :: subset
    character(128)                          :: err_message

    err_message = ""
    if(.not. is_allocated(full)) then
      err_message = "source_func_lw%subset: Asking for a subset of unallocated data"
      return
    end if
    if(start < 1 .or. start + n-1 > get_ncol(full)) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    !
    ! Could check to see if subset is correctly sized, has consistent spectral discretization
    !
    if(is_allocated(subset)) call finalize(subset)
    err_message = alloc(subset, n, get_nlay(full), full)
    if(err_message /= "") return
    subset%sfc_source    (1:n,  :) = full%sfc_source    (start:start+n-1,  :)
    subset%lay_source    (1:n,:,:) = full%lay_source    (start:start+n-1,:,:)
    subset%lev_source_inc(1:n,:,:) = full%lev_source_inc(start:start+n-1,:,:)
    subset%lev_source_dec(1:n,:,:) = full%lev_source_dec(start:start+n-1,:,:)
  end function get_subset_range_lw
  ! ------------------------------------------------------------------------------------------
  function get_subset_range_sw(full, start, n, subset) result(err_message)
    class(ty_source_func_sw), intent(inout) :: full
    integer,                  intent(in   ) :: start, n
    class(ty_source_func_sw), intent(inout) :: subset
    character(128)                          :: err_message

    err_message = ""
    if(.not. is_allocated(full)) then
      err_message = "source_func_sw%subset: Asking for a subset of unallocated data"
      return
    end if
    if(start < 1 .or. start + n-1 > get_ncol(full)) &
       err_message = "optical_props%subset: Asking for columns outside range"
    if(err_message /= "") return

    !
    ! Could check to see if subset is correctly sized, has consistent spectral discretization
    !
    if(is_allocated(subset)) call finalize(subset)
    ! Seems like I should be able to call "alloc" generically but the compilers are complaining
    err_message = alloc(subset, n, full)

    subset%toa_source(1:n,  :) = full%toa_source(start:start+n-1,  :)
  end function get_subset_range_sw
end module mo_source_functions
