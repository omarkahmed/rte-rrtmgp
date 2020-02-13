! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
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
! Encapsulates a collection of volume mixing ratios (concentrations) of gases.
!   Each concentration is associated with a name, normally the chemical formula.
!
! Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
!   (nlay and ncol are determined from the input arrays; self-consistency is enforced)
!   example:
!   error_msg = gas_concs%set_vmr('h2o', values(:,:))
!   error_msg = gas_concs%set_vmr('o3' , values(:)  )
!   error_msg = gas_concs%set_vmr('co2', value      )
!
! Values can be requested as profiles (valid only if there are no 2D fields present in the object)
!   or as 2D fields. Values for all columns are returned although the entire collection
!   can be subsetted in the column dimension
!
! Subsets can be extracted in the column dimension
!
! Functions return strings. Non-empty strings indicate an error.
!
! -------------------------------------------------------------------------------------------------

module mo_gas_concentrations
  use mo_rte_kind,           only: wp
  use mo_rrtmgp_util_string, only: lower_case, char_f2c, char_c2f
  use mo_rte_util_array,     only: any_vals_outside
  implicit none
  integer, parameter :: GAS_NOT_IN_LIST = -1

  type :: conc_field
    real(wp), dimension(:,:), pointer :: conc => NULL()
  end type conc_field

  type, public :: ty_gas_concs
    ! Data
    character(len=32), dimension(:), allocatable :: gas_name
    type(conc_field),  dimension(:), allocatable :: concs
    integer :: ncol = 0
    integer :: nlay = 0
  end type ty_gas_concs

  interface get_vmr
    module procedure get_vmr_1d, get_vmr_2d
  end interface

  interface set_vmr
    module procedure set_vmr_1d, set_vmr_2d, set_vmr_scalar
  end interface


contains
  ! -------------------------------------------------------------------------------------
  function init(cls, gas_names) result(error_msg)
    type(ty_gas_concs),            intent(inout) :: cls
    character(len=*), dimension(:), intent(in   ) :: gas_names
    character(len=128)                            :: error_msg
    character(len=128) :: tmpstri, tmpstrj
    ! ---------
    integer :: i, j, ngas
    ! ---------
    error_msg = ''
    ngas = size(gas_names)
    !
    ! Check for no duplicate gas names, no empty names
    !
    if(any(len_trim(gas_names) == 0)) &
      error_msg = "ty_gas_concs%init(): must provide non-empty gas names"

    do i = 1, ngas-1
      do j = i+1, ngas
        call char_f2c( gas_names(i) , tmpstri )
        call char_f2c( gas_names(j) , tmpstrj )
        call lower_case( tmpstri , tmpstri )
        call lower_case( tmpstrj , tmpstrj )
        call char_c2f( tmpstri , tmpstri )
        call char_c2f( tmpstrj , tmpstrj )
        if ( trim(tmpstri) == trim(tmpstrj) ) then
          error_msg = "ty_gas_concs%init(): duplicate gas names aren't allowed"
          exit
        end if
      end do
    end do
    if(error_msg /= "") return
    !
    ! Allocate fixed-size arrays
    !
    call reset(cls)
    allocate(cls%gas_name(ngas), cls%concs(ngas))

    cls%gas_name(:) = gas_names(:)
  end function
  ! -------------------------------------------------------------------------------------
  !
  ! Set concentrations --- scalar, 1D, 2D
  !
  ! -------------------------------------------------------------------------------------
  function set_vmr_scalar(cls, gas, w) result(error_msg)
    ! In OpenACC context scalar w always assumed to be on the CPU
    type(ty_gas_concs), intent(inout) :: cls
    character(len=*),    intent(in   ) :: gas
    real(wp),            intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''
    if (w < 0._wp .or. w > 1._wp) then
      error_msg = 'ty_gas_concs%set_vmr(): concentrations should be >= 0, <= 1'
      return
    endif

    igas = find_gas(cls,gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // ' but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(cls%concs(igas)%conc)) then
      if ( any(shape(cls%concs(igas)%conc) /= [1, 1]) ) then
        deallocate(cls%concs(igas)%conc)
        nullify   (cls%concs(igas)%conc)
      end if
    end if
    if (.not. associated(cls%concs(igas)%conc)) then
      allocate(cls%concs(igas)%conc(1,1))
    end if

    cls%concs(igas)%conc(:,:) = w
  end function set_vmr_scalar
  ! -------------------------------------------------------------------------------------
  function set_vmr_1d(cls, gas, w) result(error_msg)
    ! In OpenACC context w assumed to be either on the CPU or on the GPU
    type(ty_gas_concs), intent(inout) :: cls
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:), &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any_vals_outside(w, 0._wp, 1._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1'
    endif
    if(cls%nlay > 0) then
      if(size(w) /= cls%nlay) error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      cls%nlay = size(w)
    end if
    if(error_msg /= "") return

    igas = find_gas(cls,gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // ' but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(cls%concs(igas)%conc)) then
      if ( any(shape(cls%concs(igas)%conc) /= [1, cls%nlay]) ) then
        deallocate(cls%concs(igas)%conc)
        nullify   (cls%concs(igas)%conc)
      end if
    end if
    if (.not. associated(cls%concs(igas)%conc)) then
      allocate(cls%concs(igas)%conc(1,cls%nlay))
    end if

    cls%concs(igas)%conc(1,:) = w

  end function set_vmr_1d
  ! -------------------------------------------------------------------------------------
  function set_vmr_2d(cls, gas, w) result(error_msg)
    ! In OpenACC context w assumed to be either on the CPU or on the GPU
    type(ty_gas_concs), intent(inout) :: cls
    character(len=*),    intent(in   ) :: gas
    real(wp), dimension(:,:),  &
                         intent(in   ) :: w
    character(len=128)                 :: error_msg
    ! ---------
    integer :: igas
    ! ---------
    error_msg = ''

    if (any_vals_outside(w, 0._wp, 1._wp)) then
      error_msg = 'ty_gas_concs%set_vmr: concentrations should be >= 0, <= 1'
    endif

    if(cls%ncol > 0 .and. size(w, 1) /= cls%ncol) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (ncol)'
    else
      cls%ncol = size(w, 1)
    end if

    if(cls%nlay > 0 .and. size(w, 2) /= cls%nlay) then
      error_msg = 'ty_gas_concs%set_vmr: different dimension (nlay)'
    else
      cls%nlay = size(w, 2)
    end if
    if(error_msg /= "") return

    igas = find_gas(cls,gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%set_vmr(): trying to set ' // trim(gas) // 'but name not provided at initialization'
      return
    end if
    !
    ! Deallocate anything existing -- could be more efficient to test if it's already the correct size
    !
    ! This cannot be made a function, because we need all the hierarchy for the correct OpenACC attach
    if (associated(cls%concs(igas)%conc)) then
      if ( any(shape(cls%concs(igas)%conc) /= [cls%ncol,cls%nlay]) ) then
        deallocate(cls%concs(igas)%conc)
        nullify   (cls%concs(igas)%conc)
      end if
    end if
    if (.not. associated(cls%concs(igas)%conc)) then
      allocate(cls%concs(igas)%conc(cls%ncol,cls%nlay))
    end if

    cls%concs(igas)%conc(:,:) = w(:,:)
  end function set_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Return volume mixing ratio as 1D or 2D array
  !
  ! -------------------------------------------------------------------------------------
  !
  ! 1D array ( lay depdendence only)
  !
  function get_vmr_1d(cls, gas, array) result(error_msg)
    type(ty_gas_concs) :: cls
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:),   intent(out) :: array
    character(len=128) :: error_msg
    ! ---------------------
    integer :: igas
    ! ---------------------
    error_msg = ''

    igas = find_gas(cls,gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
    else if(.not. associated(cls%concs(igas)%conc)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // " concentration hasn't been set"
    else if(size(cls%concs(igas)%conc, 1) > 1) then ! Are we requesting a single profile when many are present?
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' requesting single profile but many are available'
    end if

    if(cls%nlay > 0 .and. cls%nlay /= size(array)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
    end if
    if(error_msg /= "") return

    if(size(cls%concs(igas)%conc, 2) > 1) then
      array(:) = cls%concs(igas)%conc(1,:)
    else
      array(:) = cls%concs(igas)%conc(1,1)
    end if

  end function get_vmr_1d
  ! -------------------------------------------------------------------------------------
  !
  ! 2D array (col, lay)
  !
  function get_vmr_2d(cls, gas, array) result(error_msg)
    type(ty_gas_concs) :: cls
    character(len=*),         intent(in ) :: gas
    real(wp), dimension(:,:), intent(out) :: array
    character(len=128)                    :: error_msg
    ! ---------------------
    integer :: icol, ilay, igas
    ! ---------------------
    error_msg = ''

    igas = find_gas(cls,gas)
    if (igas == GAS_NOT_IN_LIST) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' not found'
    else if(.not. associated(cls%concs(igas)%conc)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // " concentration hasn't been set"
    end if
    !
    ! Is the requested array the correct size?
    !
    if(cls%ncol > 0 .and. cls%ncol /= size(array,1)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (ncol)'
    end if
    if(cls%nlay > 0 .and. cls%nlay /= size(array,2)) then
      error_msg = 'ty_gas_concs%get_vmr; gas ' // trim(gas) // ' array is wrong size (nlay)'
    end if
    if(error_msg /= "") return

    if(size(cls%concs(igas)%conc, 1) > 1) then      ! Concentration stored as 2D
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
          !print *, (size(cls%concs))
          array(icol,ilay) = cls%concs(igas)%conc(icol,ilay)
        end do
      end do
    else if(size(cls%concs(igas)%conc, 2) > 1) then ! Concentration stored as 1D
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
         array(icol, ilay) = cls%concs(igas)%conc(1,ilay)
        end do
      end do
    else                                             ! Concentration stored as scalar
      do ilay = 1, size(array,2)
        do icol = 1, size(array,1)
          array(icol,ilay) = cls%concs(igas)%conc(1,1)
        end do
      end do
    end if

  end function get_vmr_2d
  ! -------------------------------------------------------------------------------------
  !
  ! Extract a subset of n columns starting with column 'start'
  !
  ! -------------------------------------------------------------------------------------
  function get_subset(cls, start, n, subset) result(error_msg)
    type(ty_gas_concs),      intent(in   ) :: cls
    integer,                  intent(in   ) :: start, n
    type(ty_gas_concs),      intent(inout) :: subset
    character(len=128)                      :: error_msg
    ! ---------------------
    integer :: i
    ! ---------------------
    error_msg = ''
    if(n <= 0) &
       error_msg = "gas_concs%get_vmr: Asking for 0 or fewer columns "
    if(start < 1 ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(cls%ncol > 0 .and. start > cls%ncol .or. start+n-1 > cls%ncol ) &
       error_msg = "gas_concs%get_vmr: Asking for columns outside range"
    if(error_msg /= "") return

    call reset(subset)
    allocate(subset%gas_name(size(cls%gas_name)), &
             subset%concs   (size(cls%concs))) ! These two arrays should be the same length
    subset%nlay = cls%nlay
    subset%ncol = merge(n, 0, cls%ncol > 0)
    subset%gas_name(:)  = cls%gas_name(:)

    do i = 1, size(cls%gas_name)
      !
      ! Preserve scalar/1D/2D representation in subset,
      !   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits cls dimension)
      !
      allocate(subset%concs(i)%conc(min(max(subset%ncol,1), size(cls%concs(i)%conc, 1)), &
                                    min(    subset%nlay,    size(cls%concs(i)%conc, 2))))
      if(size(cls%concs(i)%conc, 1) > 1) then      ! Concentration stored as 2D
        subset%concs(i)%conc(:,:) = cls%concs(i)%conc(start:(start+n-1),:)
      else
        subset%concs(i)%conc(:,:) = cls%concs(i)%conc(:,:)
      end if
    end do

  end function get_subset
  ! -------------------------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  ! -------------------------------------------------------------------------------------
  subroutine reset(cls)
    type(ty_gas_concs), intent(inout) :: cls
    ! -----------------
    integer :: i
    ! -----------------
    cls%nlay = 0
    cls%ncol = 0
    if(allocated(cls%gas_name)) deallocate(cls%gas_name)
    if (allocated(cls%concs)) then
      do i = 1, size(cls%concs)
        if(associated(cls%concs(i)%conc)) then
          deallocate(cls%concs(i)%conc)
          nullify(cls%concs(i)%conc)
        end if
      end do
      deallocate(cls%concs)
    end if
  end subroutine reset
  ! -------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  ! -------------------------------------------------------------------------------------
  function get_num_gases(cls)
    type(ty_gas_concs), intent(in) :: cls
    integer :: get_num_gases

    get_num_gases = size(cls%gas_name)
    return
  end function get_num_gases
  ! -------------------------------------------------------------------------------------
  function get_gas_names(cls)
    type(ty_gas_concs), intent(in) :: cls
    character(len=32), dimension(size(cls%gas_name)) :: get_gas_names

    get_gas_names(:) = cls%gas_name(:)
    return
  end function get_gas_names
  ! -------------------------------------------------------------------------------------
  !
  ! Private procedures
  !
  ! -------------------------------------------------------------------------------------
  !
  ! find gas in list; GAS_NOT_IN_LIST if not found
  !
  function find_gas(cls, gas)
    character(len=*),    intent(in) :: gas
    type(ty_gas_concs), intent(in) :: cls
    integer                         :: find_gas
    ! -----------------
    integer :: igas
    character(len=128) :: tmpstr1, tmpstr2
    ! -----------------
    find_gas = GAS_NOT_IN_LIST
    if(.not. allocated(cls%gas_name)) return
    do igas = 1, size(cls%gas_name)
      call char_f2c( cls%gas_name(igas) , tmpstr1 )
      call char_f2c( gas                 , tmpstr2 )
      call lower_case( tmpstr1 , tmpstr1 )
      call lower_case( tmpstr2 , tmpstr2 )
      call char_c2f( tmpstr1 , tmpstr1 )
      call char_c2f( tmpstr2 , tmpstr2 )
      if ( trim(tmpstr1) == trim(tmpstr2) ) then
        find_gas = igas
      end if
    end do
  end function
end module mo_gas_concentrations
