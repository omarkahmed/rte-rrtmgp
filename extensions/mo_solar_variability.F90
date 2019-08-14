!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2016-2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Optional calculation of solar variability facular and 
! sunspot indices

module mo_solar_variability
  use mo_rte_kind,           only: wp
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp

!  use mo_rrtmgp_constants, only: cp_dry, grav ! Only needed for heating rate calculation
  implicit none
  private
  type, public :: ty_solar_var
      !
      ! Data
      ! 
      real(wp), dimension(:,:), allocatable :: avgcyc_ind  ! solar variabilty index lookup table 
                                                           ! time-averaged over solar cycles 13-24.
                                                           ! (NRLSSI2 facular "Bremen" index and 
                                                           ! sunspot "SPOT67" index)
                                                           ! (nsolarterms, nsolarfrac) -> (2,134)
    contains
      !
      ! Public procedures
      !
      generic,   public  :: solar_var_index => solar_var_ind_direct, &
                                               solar_var_ind_interp 
      procedure, public  :: set_scon_or_tsi
      procedure, public  :: load_avgcyc
      procedure, public  :: finalize
      !
      ! Internal procedures
      !
      procedure, private :: solar_var_ind_direct
      procedure, private :: solar_var_ind_interp
      !
      ! Public interface
      !
  end type ty_solar_var

  interface check_range
    module procedure check_range_scalar
  end interface check_range

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routine to load mean facular and sunspot index tables
  !
  ! ------------------------------------------------------------------------------
  function set_scon_or_tsi(this, gas_opt, scon_or_tsi) result(error_msg)

    class(ty_solar_var),        intent(in   ) :: this
    type(ty_gas_optics_rrtmgp), intent(inout) :: gas_opt     ! Gas optics type
    ! Lookup table of mean solar cycle facular brightening and sunspot dimming indices
    real(wp), optional,         intent(in   ) :: scon_or_tsi
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    real(wp)              :: scon_def

    ! Default solar constant or TSI (1360.858 Wm-2)
    if (allocated(gas_opt%solar_irr_int)) then 
       scon_def = gas_opt%solar_irr_int(1) + gas_opt%solar_irr_int(2) + &
                  gas_opt%solar_irr_int(3)
       gas_opt%scon_or_tsi = scon_def
    endif
    ! User-specified solar constant or TSI if different from default
    if (present(scon_or_tsi) .and.  scon_or_tsi /= scon_def) then 
       gas_opt%scon_or_tsi = scon_or_tsi
    endif

    error_msg = check_range(gas_opt%scon_or_tsi, 0._wp, huge(gas_opt%scon_or_tsi), &
                            'set_scon_or_tsi: scon_or_tsi must be positive')
  
  end function set_scon_or_tsi
  ! ------------------------------------------------------------------------------
  !
  ! Routine to load mean facular and sunspot index tables
  !
  ! ------------------------------------------------------------------------------
  function load_avgcyc(this, avgcyc_ind) result(error_msg)

    class(ty_solar_var),      intent(inout) :: this
    ! Lookup table of mean solar cycle facular brightening and sunspot dimming indices
    real(wp), dimension(:,:), intent(in   ) :: avgcyc_ind   
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer               :: nsolarterms, nsolarfrac

    error_msg = ""
    !
    ! LUT index dimensions
    !
    nsolarterms= size(avgcyc_ind,dim=1)
    nsolarfrac = size(avgcyc_ind,dim=2)
    !
    ! Allocate LUT index array
    allocate(this%avgcyc_ind(nsolarterms, nsolarfrac))

    ! Load LUT index array
    this%avgcyc_ind = avgcyc_ind

  end function load_avgcyc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Finalize
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize(this)
    class(ty_solar_var), intent(inout) :: this

    ! Lookup table solar variability indices
    if(allocated(this%avgcyc_ind)) then
      deallocate(this%avgcyc_ind)
    end if

  end subroutine finalize
  !------------------------------------------------------------------------------------------
  !
  ! Provides the user two options for setting the facular brightening and sunspot dimming 
  ! indices required to allow solar variability. Only one option should be selected. 
  ! 
  ! Option 1: solar_var_ind_direct
  ! Direct specification of the facular and sunspot indices by the user. 
  !
  function solar_var_ind_direct(this,                      &
                                mg_index, sb_index,        &
                                gas_opt)                   &
                                result(error_msg)

    class(ty_solar_var),        intent(in   ) :: this        ! Solar variability
    real(wp),                   intent(in   ) :: mg_index, & ! Facular brightening index (user-defined)
                                                 sb_index    ! Sunspot dimmng index (user-defined)
    type(ty_gas_optics_rrtmgp), intent(inout) :: gas_opt     ! Gas optics type

    character(len=128)                        :: error_msg

    ! ----------------------------------------------------------
    ! Local variables
    !
    ! ----------------------------------------------------------
    !
    ! Error checking

    ! Check input data values
    !
    error_msg = check_range(mg_index, 0._wp, huge(mg_index), 'mg_index')
    if(error_msg  /= '') return
    error_msg = check_range(sb_index, 0._wp, huge(sb_index), 'sb_index')
    if(error_msg  /= '') return
    !
    ! Allocate solar variabilty index array
    allocate(gas_opt%solar_var_ind(2))
    !
    ! Direct specification of facular and sunspot indices
    !
    gas_opt%solar_var_ind(1) = mg_index
    gas_opt%solar_var_ind(2) = sb_index

  end function solar_var_ind_direct
  ! 
  ! Option 2: solar_var_ind_interp
  ! Facular brightening and sunspot dimming indices are derived from the 
  ! averaged solar cycle, which is the mean of Solar Cycles 13-24. The user specifices 
  ! the solar cycle fraction (0 to 1) and the indices are interpolated to the
  ! requested fractional position within the cycle, where 0 is close to solar minimum.
  ! 
  function solar_var_ind_interp(this,                   &
                                solcycfrac,             &
                                gas_opt,                &
                                mg_scl, sb_scl)         &
                                result(error_msg)

    class(ty_solar_var),        intent(in   ) :: this        ! Solar variability
    real(wp),                   intent(in   ) :: solcycfrac  ! solar cycle fraction
    type(ty_gas_optics_rrtmgp), intent(inout) :: gas_opt     ! Gas optics type

    real(wp), optional,         intent(in   ) :: mg_scl  , & ! Facular brightening scaling (user-defined)
                                                 sb_scl      ! Sunspot dimmng scaling (user-defined)

    character(len=128)                        :: error_msg

    ! ----------------------------------------------------------
    ! Local variables
    !
    integer  :: nsolfrac                                 ! Number of solar fraction points in facular 
                                                         ! and sunspot tables
    integer  :: sfid                                     ! Solar variability solar cycle fraction index

    real(wp) :: svar_fac_0, svar_spt_0                   ! Interpolated indices
    real(wp) :: intrvl_len                               ! Fractional interval length of mgavgcyc and sbavgcyc
    real(wp) :: intrvl_len_hf                            ! Fractional half interval length of mgavgcyc and sbavgcyc
    real(wp) :: fraclo, frachi, intfrac                  ! Solar variability interpolation factors

    real(wp), parameter ::  solcycfrac_min = 0.0189_wp   ! Solar cycle fraction at mean solar cycle minimum
    real(wp), parameter ::  solcycfrac_max = 0.3750_wp   ! Solar cycle fraction at mean solar cycle maximum
    real(wp), parameter ::  fracdiff_min2max = 0.3561_wp ! 0.3750 - 0.0189
    real(wp), parameter ::  fracdiff_max2min = 0.6439_wp ! 1.0189 - 0.3750
    real(wp) :: wgt                                      ! Weighting factor for amplitude scale factor adjustment
    ! ----------------------------------------------------------
    !
    ! Error checking
    !
    ! Check input data sizes and values
    !
    error_msg = check_range(solcycfrac, 0._wp, 1._wp, 'solcycfrac')
    if(error_msg  /= '') return
    !
    ! Interpolate solar variability indices to requested solar cycle fraction,
    ! and derive final facular and sunspot indices, with optional scaling
    !
    nsolfrac = size(this%avgcyc_ind,2)
    if (solcycfrac .le. 0._wp) then
       svar_fac_0 = this%avgcyc_ind(1,1)
       svar_spt_0 = this%avgcyc_ind(2,1)
    elseif (solcycfrac .ge. 1._wp) then
       svar_fac_0 = this%avgcyc_ind(1,nsolfrac)
       svar_spt_0 = this%avgcyc_ind(2,nsolfrac)
    else
       intrvl_len = 1._wp / (nsolfrac-2)
       intrvl_len_hf = 0.5_wp * intrvl_len
       !   Initial half interval (1)
       if (solcycfrac .le. intrvl_len_hf) then 
          sfid = 1
          fraclo = 0._wp
          frachi = intrvl_len_hf
       endif
       !   Main whole intervals (131)
       if (solcycfrac .gt. intrvl_len_hf .and. solcycfrac .lt. 1._wp-intrvl_len_hf) then 
          sfid = floor((solcycfrac-intrvl_len_hf) * (nsolfrac-2)) + 2
          fraclo = (sfid-2) * intrvl_len + intrvl_len_hf
          frachi = fraclo + intrvl_len
       endif
       !   Final half interval (1)
       if (solcycfrac .ge. 1._wp-intrvl_len_hf) then 
          sfid = (nsolfrac-2) + 1
          fraclo = 1._wp - intrvl_len_hf
          frachi = 1._wp
       endif
       intfrac = (solcycfrac - fraclo) / (frachi - fraclo)
       svar_fac_0 = this%avgcyc_ind(1,sfid) + &
                    intfrac * (this%avgcyc_ind(1,sfid+1) - this%avgcyc_ind(1,sfid))
       svar_spt_0 = this%avgcyc_ind(2,sfid) + &
                    intfrac * (this%avgcyc_ind(2,sfid+1) - this%avgcyc_ind(2,sfid))
    endif
    !
    ! Allocate solar variabilty index array
    allocate(gas_opt%solar_var_ind(2))
    allocate(gas_opt%solar_var_ind_scl(2))
    !
    ! Adjust amplitude scaling of mean solar cycle to be 1.0 at solar minimum (solcycfrac_min=0.0189),
    ! to be the user-provided scaling at solar maximum (solcycfrac_max=0.3750), and to vary between 
    ! those values at intervening values of solcycfrac. 
    !
    gas_opt%solar_var_ind_scl(1:2) = 1._wp
    if (present(mg_scl) .and. present(sb_scl)) then 
       if (mg_scl .ne. 1._wp .or. sb_scl .ne. 1._wp) then 
          if (solcycfrac .ge. 0._wp .and. solcycfrac .lt. solcycfrac_min) then
             wgt = (solcycfrac+1._wp-solcycfrac_max)/fracdiff_max2min
             gas_opt%solar_var_ind_scl(1) = mg_scl + wgt * (1._wp-mg_scl)
             gas_opt%solar_var_ind_scl(2) = sb_scl + wgt * (1._wp-sb_scl)
          endif
          if (solcycfrac .ge. solcycfrac_min .and. solcycfrac .le. solcycfrac_max) then
             wgt = (solcycfrac-solcycfrac_min)/fracdiff_min2max
             gas_opt%solar_var_ind_scl(1) = 1._wp + wgt * (mg_scl-1._wp)
             gas_opt%solar_var_ind_scl(2) = 1._wp + wgt * (sb_scl-1._wp)
          endif
          if (solcycfrac .gt. solcycfrac_max .and. solcycfrac .le. 1._wp) then
             wgt = (solcycfrac-solcycfrac_max)/fracdiff_max2min
             gas_opt%solar_var_ind_scl(1) = mg_scl + wgt * (1._wp-mg_scl)
             gas_opt%solar_var_ind_scl(2) = sb_scl + wgt * (1._wp-sb_scl)
          endif
       endif
    endif
    !
    ! Compute final facular and sunspot indices, with optional scaling
    !
    if (present(mg_scl) .and. present(sb_scl)) then 
       gas_opt%solar_var_ind(1) = gas_opt%solar_var_ind_scl(1) * &
                                  (svar_fac_0 - gas_opt%svar_offset(1)) + gas_opt%svar_offset(1)
       gas_opt%solar_var_ind(2) = gas_opt%solar_var_ind_scl(2) * &
                                  (svar_spt_0 - gas_opt%svar_offset(2)) + gas_opt%svar_offset(2)
    else
       gas_opt%solar_var_ind(1) = (svar_fac_0 - gas_opt%svar_offset(1)) + gas_opt%svar_offset(1)
       gas_opt%solar_var_ind(2) = (svar_spt_0 - gas_opt%svar_offset(2)) + gas_opt%svar_offset(2)
    end if
 
    ! Error checking
    ! Check interpolated index values
    !
    error_msg = check_range(gas_opt%solar_var_ind(1), 0._wp, huge(gas_opt%solar_var_ind(1)), 'mg_index')
    if(error_msg  /= '') return
    error_msg = check_range(gas_opt%solar_var_ind(2), 0._wp, huge(gas_opt%solar_var_ind(2)), 'sb_index')
    if(error_msg  /= '') return

  end function solar_var_ind_interp

  ! --------------------------------------------------------------------------------------
  !
  ! Values
  !
  ! --------------------------------------------------------------------------------------
  function check_range_scalar(val, minV, maxV, label)
    real(wp),                   intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label
    character(len=128)                     :: check_range_scalar

    check_range_scalar = ""
    if(val < minV .or. val > maxV) &
      check_range_scalar = trim(label) // ' value out of range.'
  end function check_range_scalar
  ! --------------------------------------------------------------------------------------
end module mo_solar_variability
