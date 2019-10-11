! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Description: Optional calculation of solar variability facular and
! sunspot indices

module mo_solar_variability
  use mo_rte_kind,           only: wp

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
      procedure, public :: solar_var_ind_interp
      procedure, public :: load_avgcyc
      procedure, public :: finalize
      !
  end type ty_solar_var

contains
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
  ! ------------------------------------------------------------------------------
  !
  ! Facular brightening and sunspot dimming indices are derived from the
  ! averaged solar cycle, which is the mean of Solar Cycles 13-24. The user specifices
  ! the solar cycle fraction (0 to 1) and the indices are interpolated to the
  ! requested fractional position within the cycle, where 0 is close to solar minimum.
  !
  ! Optional amplitude scaling values are adjusted to the requested solar cycle
  ! fraction such that the full scaling is applied at solar maximum, no scaling
  ! is applied at solar minimum, and intermediate values are applied at other points
  ! in the solar cycle. In other workds, the amplitude of the solar cycle is stretched
  ! to match the requested scaling at solar maximum.
  !
  function solar_var_ind_interp(this,                      &
                                solcycfrac, scon,          &
                                irr_int, fac_int, spt_int, &
                                fac_avg_ind, spt_avg_ind,  &
                                fac_offset, spt_offset,    &
                                mg_index, sb_index,        &
                                tsi)                       &
                                result(error_msg)

    class(ty_solar_var),        intent(in   ) :: this        ! Solar variability
    real(wp),                   intent(in   ) :: solcycfrac  ! solar cycle fraction
    real(wp), optional,         intent(in   ) :: scon        ! Solar constant (Wm-2)
    real(wp),                   intent(in   ) :: irr_int     ! integrated quiet sun irradiance term
    real(wp),                   intent(in   ) :: fac_int     ! integrated facular brightening irradiance
    real(wp),                   intent(in   ) :: spt_int     ! integrated sunspot dimming irradiance
    real(wp),                   intent(in   ) :: fac_avg_ind ! facular brightening average index for mean solar cycle
    real(wp),                   intent(in   ) :: spt_avg_ind ! sunspot dimming average index for mean solar cycle
    real(wp),                   intent(in   ) :: fac_offset  ! facular brightening offset term
    real(wp),                   intent(in   ) :: spt_offset  ! sunspot dimming offset term
    real(wp),                   intent(out  ) :: mg_index    ! Facular brightening NRLSSI2 index
                                                             ! interpolated from the mean solar cycle
                                                             ! to the provided solar cycle fraction
    real(wp),                   intent(out  ) :: sb_index    ! Sunspot dimmng NRLSSI2 index
                                                             ! interpolated from the mean solar cycle
                                                             ! to the provided solar cycle fraction
    real(wp),                   intent(out  ) :: tsi         ! Total solar irradiance (in Wm-2) for 
                                                             ! either the default solar constant or the 
                                                             ! optional, user-provided solar constant 
                                                             ! adjusted to the derived facular and sunspot 
                                                             ! index values.

    character(len=128)                        :: error_msg

    ! ----------------------------------------------------------
    ! Local variables
    !
    integer  :: nsolfrac                                 ! Number of solar fraction points in facular
                                                         ! and sunspot tables
    integer  :: sfid                                     ! Solar variability solar cycle fraction index

    real(wp) :: intrvl_len                               ! Fractional interval length of mgavgcyc and sbavgcyc
    real(wp) :: intrvl_len_hf                            ! Fractional half interval length of mgavgcyc and sbavgcyc
    real(wp) :: fraclo, frachi, intfrac                  ! Solar variability interpolation factors
    real(wp) :: a_fac, b_spt, c_irr                         
    ! ----------------------------------------------------------
    !
    ! Error checking
    error_msg = ""
    !
    ! Check input data sizes and values
    !
    if (solcycfrac .lt. 0._wp .or. solcycfrac .gt. 1._wp) &
       error_msg = 'solar_var_ind_interp: solcycfrac out of range'
    if(error_msg  /= '') return
    !
    ! Interpolate solar variability indices to requested solar cycle fraction,
    ! and derive final facular and sunspot indices
    !
    ! nsolfrac is the length of the time dimension of the interpolation tables 
    ! of facular and sunspot indices over the mean solar cycle (this%avgcyc_ind).
    ! The end-points of avgcyc_ind represent the indices at solcycfrac values of 
    ! 0 (first day of the first year) and 1 (last day of the 11th year), while 
    ! the intervening values of avgcyc_ind represent the indices at the center 
    ! of each month over the mean 11-year solar cycle. 
    if (allocated (this%avgcyc_ind)) then
       nsolfrac = size(this%avgcyc_ind,2)
    ! Define indices for the lowest allowable value of solcycfrac
       if (solcycfrac .eq. 0._wp) then
          mg_index = this%avgcyc_ind(1,1)
          sb_index = this%avgcyc_ind(2,1)
    ! Define indices for the highest allowable value of solcycfrac
       elseif (solcycfrac .eq. 1._wp) then
          mg_index = this%avgcyc_ind(1,nsolfrac)
          sb_index = this%avgcyc_ind(2,nsolfrac)
    ! Define indices for intervening values of solcycfrac
       else
          intrvl_len = 1._wp / (nsolfrac-2)
          intrvl_len_hf = 0.5_wp * intrvl_len
    ! Define interpolation fractions for the first interval, which represents 
    ! the first half of the first month of the first year of the mean 11-year
    ! solar cycle
          if (solcycfrac .le. intrvl_len_hf) then
             sfid = 1
             fraclo = 0._wp
             frachi = intrvl_len_hf
          endif
    ! Define interpolation fractions for the intervening intervals, which represent 
    ! the center point of each month in each year of the mean 11-year solar cycle
          if (solcycfrac .gt. intrvl_len_hf .and. solcycfrac .lt. 1._wp-intrvl_len_hf) then
             sfid = floor((solcycfrac-intrvl_len_hf) * (nsolfrac-2)) + 2
             fraclo = (sfid-2) * intrvl_len + intrvl_len_hf
             frachi = fraclo + intrvl_len
          endif
    ! Define interpolation fractions for the last interval, which represents 
    ! the last half of the last month of the last year of the mean 11-year
    ! solar cycle
          if (solcycfrac .ge. 1._wp-intrvl_len_hf) then
             sfid = (nsolfrac-2) + 1
             fraclo = 1._wp - intrvl_len_hf
             frachi = 1._wp
          endif
    ! Interpolate the facular (mg_index) and sunspot (sb_index) indices for the
    ! requested value of solcycfrac
          intfrac = (solcycfrac - fraclo) / (frachi - fraclo)
          mg_index = this%avgcyc_ind(1,sfid) + &
                     intfrac * (this%avgcyc_ind(1,sfid+1) - this%avgcyc_ind(1,sfid))
          sb_index = this%avgcyc_ind(2,sfid) + &
                     intfrac * (this%avgcyc_ind(2,sfid+1) - this%avgcyc_ind(2,sfid))
       endif
    endif

  ! Derive tsi for derived facular and sunspot indices if scon provided
    a_fac = (mg_index - fac_offset) / (fac_avg_ind - fac_offset)
    b_spt = (sb_index - spt_offset) / (spt_avg_ind - spt_offset)
    c_irr = 1.0_wp
    if (present(scon)) then
       c_irr = (scon - fac_int - spt_int) / irr_int
    endif
    tsi = c_irr * irr_int + a_fac * fac_int + b_spt * spt_int

  end function solar_var_ind_interp
  ! --------------------------------------------------------------------------------------
end module mo_solar_variability
