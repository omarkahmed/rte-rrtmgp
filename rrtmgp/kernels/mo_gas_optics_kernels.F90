! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths,
!   source functions.

module mo_gas_optics_kernels
  use mo_rte_kind,      only : wp, wl
  implicit none

  interface

    subroutine interpolation( &
                  ncol,nlay,ngas,nflav,neta, npres, ntemp, &
                  flavor,                                  &
                  press_ref_log, temp_ref,press_ref_log_delta,    &
                  temp_ref_min,temp_ref_delta,press_ref_trop_log, &
                  vmr_ref,                                        &
                  play,tlay,col_gas,                              &
                  jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress) bind(C, name="interpolation")
      use mo_rte_kind,      only : wp, wl
      implicit none
      integer, value,                     intent(in) :: ncol,nlay
      integer, value,                     intent(in) :: ngas,nflav,neta,npres,ntemp
      integer,     dimension(2,nflav),    intent(in) :: flavor
      real(wp),    dimension(npres),      intent(in) :: press_ref_log
      real(wp),    dimension(ntemp),      intent(in) :: temp_ref
      real(wp), value,                    intent(in) :: press_ref_log_delta, &
                                                        temp_ref_min, temp_ref_delta, &
                                                        press_ref_trop_log
      real(wp),    dimension(2,0:ngas,ntemp), intent(in) :: vmr_ref
      real(wp),    dimension(ncol,nlay),        intent(in) :: play, tlay
      real(wp),    dimension(ncol,nlay,0:ngas), intent(in) :: col_gas
      integer,     dimension(ncol,nlay), intent(out) :: jtemp, jpress
      logical(wl), dimension(ncol,nlay), intent(out) :: tropo
      integer,     dimension(2,    nflav,ncol,nlay), intent(out) :: jeta
      real(wp),    dimension(2,    nflav,ncol,nlay), intent(out) :: col_mix
      real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(out) :: fmajor
      real(wp),    dimension(2,2,  nflav,ncol,nlay), intent(out) :: fminor
    end subroutine

    subroutine combine_and_reorder_2str(ncol, nlay, ngpt, tau_abs, tau_rayleigh, tau, ssa, g) bind(C, name="combine_and_reorder_2str")
      use mo_rte_kind,      only : wp, wl
      implicit none
      integer, value,                      intent(in) :: ncol, nlay, ngpt
      real(wp), dimension(ngpt,nlay,ncol), intent(in   ) :: tau_abs, tau_rayleigh
      real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa, g ! inout because components are allocated
    end subroutine

    function interpolate3D(scaling, fmajor, k, igpt, jeta, jtemp, jpress,ngpt,neta,npres,ntemp) result(res) bind(C, name="interpolate3D")
      use mo_rte_kind,      only : wp, wl
      implicit none
      real(wp), dimension(2),     intent(in) :: scaling
      real(wp), dimension(2,2,2), intent(in) :: fmajor ! interpolation fractions for major species
                                                       ! index(1) : reference eta level (temperature dependent)
                                                       ! index(2) : reference pressure level
                                                       ! index(3) : reference temperature level
      real(wp), dimension(ngpt,neta,npres+1,ntemp),intent(in) :: k
      integer, value,              intent(in) :: igpt, ngpt, neta, npres, ntemp
      integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
      integer, value,              intent(in) :: jtemp ! interpolation index for temperature
      integer, value,              intent(in) :: jpress ! interpolation index for pressure
      real(wp)                                :: res ! the result
    end function

    function interpolate2D(fminor, k, igpt, jeta, jtemp, ngpt, neta, ntemp) result(res) bind(C, name="interpolate2D")
      use mo_rte_kind,      only : wp, wl
      implicit none
      real(wp), dimension(2,2), intent(in) :: fminor ! interpolation fractions for minor species
                                         ! index(1) : reference eta level (temperature dependent)
                                         ! index(2) : reference temperature level
      real(wp), dimension(ngpt,neta,ntemp), intent(in) :: k ! (g-point, eta, temp)
      integer, value,             intent(in) :: igpt, jtemp, ngpt, neta, ntemp ! interpolation index for temperature
      integer, dimension(2),      intent(in) :: jeta ! interpolation index for binary species parameter (eta)
      real(wp)                             :: res ! the result
    end function

    subroutine interpolate1D(val, offset, delta, table, res, tab_d1, tab_d2) bind(C, name="interpolate1D")
      use mo_rte_kind,      only : wp, wl
      implicit none
      real(wp), value, intent(in) :: val,    & ! axis value at which to evaluate table
                                     offset, & ! minimum of table axis
                                     delta     ! step size of table axis
      integer, value, intent(in) :: tab_d1, tab_d2
      real(wp), dimension(tab_d1, tab_d2), &
                intent(in) :: table ! dimensions (axis, values)
      real(wp), intent(out) ,dimension(tab_d2) :: res
    end subroutine

    subroutine compute_Planck_source(                        &
                      ncol, nlay, nbnd, ngpt,                &
                      nflav, neta, npres, ntemp, nPlanckTemp,&
                      tlay, tlev, tsfc, sfc_lay,             &
                      fmajor, jeta, tropo, jtemp, jpress,    &
                      gpoint_bands, band_lims_gpt,           &
                      pfracin, temp_ref_min, totplnk_delta, totplnk, gpoint_flavor, &
                      sfc_src, lay_src, lev_src_inc, lev_src_dec) bind(C, name="compute_Planck_source")
      use mo_rte_kind,      only : wp, wl
      implicit none
      integer, value,                             intent(in) :: ncol, nlay, nbnd, ngpt
      integer, value,                             intent(in) :: nflav, neta, npres, ntemp, nPlanckTemp
      real(wp),    dimension(ncol,nlay  ),        intent(in) :: tlay
      real(wp),    dimension(ncol,nlay+1),        intent(in) :: tlev
      real(wp),    dimension(ncol       ),        intent(in) :: tsfc
      integer, value,                             intent(in) :: sfc_lay
      real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
      integer,     dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
      logical(wl), dimension(            ncol,nlay), intent(in) :: tropo
      integer,     dimension(            ncol,nlay), intent(in) :: jtemp, jpress
      integer, dimension(ngpt),                     intent(in) :: gpoint_bands ! start and end g-point for each band
      integer, dimension(2, nbnd),                  intent(in) :: band_lims_gpt ! start and end g-point for each band
      real(wp), value,                              intent(in) :: temp_ref_min, totplnk_delta
      real(wp), dimension(ngpt,neta,npres+1,ntemp), intent(in) :: pfracin
      real(wp), dimension(nPlanckTemp,nbnd),        intent(in) :: totplnk
      integer,  dimension(2,ngpt),                  intent(in) :: gpoint_flavor
      real(wp), dimension(ngpt,     ncol), intent(out) :: sfc_src
      real(wp), dimension(ngpt,nlay,ncol), intent(out) :: lay_src
      real(wp), dimension(ngpt,nlay,ncol), intent(out) :: lev_src_inc, lev_src_dec
    end subroutine

  end interface


contains


  ! --------------------------------------------------------------------------------------
  !
  ! Compute minor and major species opitcal depth from pre-computed interpolation coefficients
  !   (jeta,jtemp,jpress)
  !
  subroutine compute_tau_absorption(                &
                ncol,nlay,nbnd,ngpt,                &  ! dimensions
                ngas,nflav,neta,npres,ntemp,        &
                nminorlower, nminorklower,          & ! number of minor contributors, total num absorption coeffs
                nminorupper, nminorkupper,          &
                idx_h2o,                            &
                gpoint_flavor,                      &
                band_lims_gpt,                      &
                kmajor,                             &
                kminor_lower,                       &
                kminor_upper,                       &
                minor_limits_gpt_lower,             &
                minor_limits_gpt_upper,             &
                minor_scales_with_density_lower,    &
                minor_scales_with_density_upper,    &
                scale_by_complement_lower,          &
                scale_by_complement_upper,          &
                idx_minor_lower,                    &
                idx_minor_upper,                    &
                idx_minor_scaling_lower,            &
                idx_minor_scaling_upper,            &
                kminor_start_lower,                 &
                kminor_start_upper,                 &
                tropo,                              &
                col_mix,fmajor,fminor,              &
                play,tlay,col_gas,                  &
                jeta,jtemp,jpress,                  &
                tau) bind(C, name="compute_tau_absorption")
    ! ---------------------
    ! input dimensions
    integer,                                intent(in) :: ncol,nlay,nbnd,ngpt
    integer,                                intent(in) :: ngas,nflav,neta,npres,ntemp
    integer,                                intent(in) :: nminorlower, nminorklower,nminorupper, nminorkupper
    integer,                                intent(in) :: idx_h2o
    ! ---------------------
    ! inputs from object
    integer,     dimension(2,ngpt),                  intent(in) :: gpoint_flavor
    integer,     dimension(2,nbnd),                  intent(in) :: band_lims_gpt
    real(wp),    dimension(ngpt,neta,npres+1,ntemp), intent(in) :: kmajor
    real(wp),    dimension(nminorklower,neta,ntemp), intent(in) :: kminor_lower
    real(wp),    dimension(nminorkupper,neta,ntemp), intent(in) :: kminor_upper
    integer,     dimension(2,nminorlower),           intent(in) :: minor_limits_gpt_lower
    integer,     dimension(2,nminorupper),           intent(in) :: minor_limits_gpt_upper
    logical(wl), dimension(  nminorlower),           intent(in) :: minor_scales_with_density_lower
    logical(wl), dimension(  nminorupper),           intent(in) :: minor_scales_with_density_upper
    logical(wl), dimension(  nminorlower),           intent(in) :: scale_by_complement_lower
    logical(wl), dimension(  nminorupper),           intent(in) :: scale_by_complement_upper
    integer,     dimension(  nminorlower),           intent(in) :: idx_minor_lower
    integer,     dimension(  nminorupper),           intent(in) :: idx_minor_upper
    integer,     dimension(  nminorlower),           intent(in) :: idx_minor_scaling_lower
    integer,     dimension(  nminorupper),           intent(in) :: idx_minor_scaling_upper
    integer,     dimension(  nminorlower),           intent(in) :: kminor_start_lower
    integer,     dimension(  nminorupper),           intent(in) :: kminor_start_upper
    logical(wl), dimension(ncol,nlay),               intent(in) :: tropo
    ! ---------------------
    ! inputs from profile or parent function
    real(wp), dimension(2,    nflav,ncol,nlay       ), intent(in) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay       ), intent(in) :: fmajor
    real(wp), dimension(2,2,  nflav,ncol,nlay       ), intent(in) :: fminor
    real(wp), dimension(            ncol,nlay       ), intent(in) :: play, tlay      ! pressure and temperature
    real(wp), dimension(            ncol,nlay,0:ngas), intent(in) :: col_gas
    integer,  dimension(2,    nflav,ncol,nlay       ), intent(in) :: jeta
    integer,  dimension(            ncol,nlay       ), intent(in) :: jtemp
    integer,  dimension(            ncol,nlay       ), intent(in) :: jpress
    ! ---------------------
    ! output - optical depth
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! ---------------------
    ! Local variables
    !
    logical(wl)                :: top_at_1
    integer, dimension(ncol,2) :: itropo_lower, itropo_upper
    integer                    :: icol, idx_tropo

    ! ----------------------------------------------------------------

    !$acc enter data create(itropo_lower, itropo_upper)
    !$acc enter data copyin(play, tlay, tropo, gpoint_flavor, jeta, jtemp, col_gas, fminor, tau)

    ! ---------------------
    ! Layer limits of upper, lower atmospheres
    ! ---------------------
    top_at_1 = play(1,1) < play(1, nlay)
    if(top_at_1) then
      !$acc parallel loop
      do icol = 1,ncol
        itropo_lower(icol,2) = nlay
        itropo_lower(icol,1) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
        itropo_upper(icol,1) = 1
        itropo_upper(icol,2) = maxloc(play(icol,:), dim=1, mask=(.not. tropo(icol,:)))
      end do
    else
      !$acc parallel loop
      do icol = 1,ncol
        itropo_lower(icol,1) = 1
        itropo_lower(icol,2) = minloc(play(icol,:), dim=1, mask=tropo(icol,:))
        itropo_upper(icol,2) = nlay
        itropo_upper(icol,1) = maxloc(play(icol,:), dim=1, mask=(.not.tropo(icol,:)))
      end do
    end if
    ! ---------------------
    ! Major Species
    ! ---------------------
    call gas_optical_depths_major(   &
          ncol,nlay,nbnd,ngpt,       & ! dimensions
          nflav,neta,npres,ntemp,    &
          gpoint_flavor,             &
          band_lims_gpt,             &
          kmajor,                    &
          col_mix,fmajor,            &
          jeta,tropo,jtemp,jpress,   &
          tau)
    ! ---------------------
    ! Minor Species - lower
    ! ---------------------
    idx_tropo = 1
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,             & ! dimensions
           ngas,nflav,ntemp,neta,      &
           nminorlower,nminorklower,   &
           idx_h2o,idx_tropo,          &
           gpoint_flavor,              &
           kminor_lower,               &
           minor_limits_gpt_lower,     &
           minor_scales_with_density_lower, &
           scale_by_complement_lower,  &
           idx_minor_lower,            &
           idx_minor_scaling_lower,    &
           kminor_start_lower,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_lower,jtemp,         &
           tau)
    ! ---------------------
    ! Minor Species - upper
    ! ---------------------
    idx_tropo = 2
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,             & ! dimensions
           ngas,nflav,ntemp,neta,      &
           nminorupper,nminorkupper,   &
           idx_h2o,idx_tropo,          &
           gpoint_flavor,              &
           kminor_upper,               &
           minor_limits_gpt_upper,     &
           minor_scales_with_density_upper, &
           scale_by_complement_upper,  &
           idx_minor_upper,            &
           idx_minor_scaling_upper,    &
           kminor_start_upper,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_upper,jtemp,         &
           tau)

    !$acc exit data delete(itropo_lower,itropo_upper)
    !$acc exit data delete(play, tlay, tropo, gpoint_flavor, jeta, jtemp, col_gas, fminor)
    !$acc exit data copyout(tau)

  end subroutine compute_tau_absorption
  ! --------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_major(ncol,nlay,nbnd,ngpt,&
                                      nflav,neta,npres,ntemp,      & ! dimensions
                                      gpoint_flavor, band_lims_gpt,   & ! inputs from object
                                      kmajor,                         &
                                      col_mix,fmajor,                 &
                                      jeta,tropo,jtemp,jpress,        & ! local input
                                      tau) bind(C, name="gas_optical_depths_major")
    ! input dimensions
    integer, intent(in) :: ncol, nlay, nbnd, ngpt, nflav,neta,npres,ntemp  ! dimensions

    ! inputs from object
    integer,  dimension(2,ngpt),  intent(in) :: gpoint_flavor
    integer,  dimension(2,nbnd),  intent(in) :: band_lims_gpt ! start and end g-point for each band
    real(wp), dimension(ngpt,neta,npres+1,ntemp), intent(in) :: kmajor

    ! inputs from profile or parent function
    real(wp),    dimension(2,    nflav,ncol,nlay), intent(in) :: col_mix
    real(wp),    dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,     dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical(wl), dimension(ncol,nlay), intent(in) :: tropo
    integer,     dimension(ncol,nlay), intent(in) :: jtemp, jpress

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp) :: tau_major ! major species optical depth
    ! local index
    integer :: icol, ilay, iflav, igpt, itropo

    ! -----------------

    ! -----------------

    ! optical depth calculation for major species
    !$acc parallel loop collapse(3)
    do ilay = 1, nlay
      do icol = 1, ncol
        ! optical depth calculation for major species
        do igpt = 1, ngpt
          ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
          itropo = merge(1,2,tropo(icol,ilay))  ! WS: moved inside innermost loop

          ! binary species parameter (eta) and col_mix depend on band flavor
          iflav = gpoint_flavor(itropo, igpt)
          tau_major = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D(col_mix(:,iflav,icol,ilay), &
                          fmajor(:,:,:,iflav,icol,ilay), kmajor, &
                          igpt, jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo,ngpt,neta,npres,ntemp)
          tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_major
        end do ! igpt
      end do
    end do ! ilay
  end subroutine gas_optical_depths_major

  ! ----------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_minor(ncol,nlay,ngpt,        &
                                      ngas,nflav,ntemp,neta, &
                                      nminor,nminork,        &
                                      idx_h2o,idx_tropo,     &
                                      gpt_flv,               &
                                      kminor,                &
                                      minor_limits_gpt,      &
                                      minor_scales_with_density,    &
                                      scale_by_complement,   &
                                      idx_minor, idx_minor_scaling, &
                                      kminor_start,        &
                                      play, tlay,          &
                                      col_gas,fminor,jeta, &
                                      layer_limits,jtemp,  &
                                      tau) bind(C, name="gas_optical_depths_minor")
    integer,                                     intent(in   ) :: ncol,nlay,ngpt
    integer,                                     intent(in   ) :: ngas,nflav
    integer,                                     intent(in   ) :: ntemp,neta,nminor,nminork
    integer,                                     intent(in   ) :: idx_h2o, idx_tropo
    integer,     dimension(2, ngpt),             intent(in   ) :: gpt_flv
    real(wp),    dimension(nminork,neta,ntemp),  intent(in   ) :: kminor
    integer,     dimension(2,nminor),            intent(in   ) :: minor_limits_gpt
    logical(wl), dimension(  nminor),            intent(in   ) :: minor_scales_with_density
    logical(wl), dimension(  nminor),            intent(in   ) :: scale_by_complement
    integer,     dimension(  nminor),            intent(in   ) :: kminor_start
    integer,     dimension(  nminor),            intent(in   ) :: idx_minor, idx_minor_scaling
    real(wp),    dimension(ncol,nlay),           intent(in   ) :: play, tlay
    real(wp),    dimension(ncol,nlay,0:ngas),    intent(in   ) :: col_gas
    real(wp),    dimension(2,2,nflav,ncol,nlay), intent(in   ) :: fminor
    integer,     dimension(2,  nflav,ncol,nlay), intent(in   ) :: jeta
    integer,     dimension(ncol, 2),             intent(in   ) :: layer_limits
    integer,     dimension(ncol,nlay),           intent(in   ) :: jtemp
    real(wp),    dimension(ngpt,nlay,ncol),      intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp), parameter :: PaTohPa = 0.01
    real(wp) :: vmr_fact, dry_fact             ! conversion from column abundance to dry vol. mixing ratio;
    real(wp) :: scaling, kminor_loc, tau_minor ! minor species absorption coefficient, optical depth
    integer  :: icol, ilay, iflav, igpt, imnr
    integer  :: gptS, gptE
    integer  :: minor_start, minor_loc, extent

    real(wp) :: myplay, mytlay, mycol_gas_h2o, mycol_gas_imnr, mycol_gas_0
    real(wp) :: myfminor(2,2)
    integer  :: myjtemp, myjeta(2), max_gpt_diff, igpt0
    ! -----------------

    extent = size(scale_by_complement,dim=1)

    ! Find the largest number of g-points per band
    max_gpt_diff = maxval( minor_limits_gpt(2,:) - minor_limits_gpt(1,:) )

    !$acc parallel loop gang vector collapse(3)
    do ilay = 1 , nlay
      do icol = 1, ncol
        do igpt0 = 0, max_gpt_diff
          !
          ! This check skips individual columns with no pressures in range
          !
          if ( layer_limits(icol,1) <= 0 .or. ilay < layer_limits(icol,1) .or. ilay > layer_limits(icol,2) ) cycle

          myplay  = play (icol,ilay)
          mytlay  = tlay (icol,ilay)
          myjtemp = jtemp(icol,ilay)
          mycol_gas_h2o = col_gas(icol,ilay,idx_h2o)
          mycol_gas_0   = col_gas(icol,ilay,0)

          do imnr = 1, extent

            scaling = col_gas(icol,ilay,idx_minor(imnr))
            if (minor_scales_with_density(imnr)) then
              !
              ! NOTE: P needed in hPa to properly handle density scaling.
              !
              scaling = scaling * (PaTohPa * myplay/mytlay)

              if(idx_minor_scaling(imnr) > 0) then  ! there is a second gas that affects this gas's absorption
                mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr))
                vmr_fact = 1._wp / mycol_gas_0
                dry_fact = 1._wp / (1._wp + mycol_gas_h2o * vmr_fact)
                ! scale by density of special gas
                if (scale_by_complement(imnr)) then ! scale by densities of all gases but the special one
                  scaling = scaling * (1._wp - mycol_gas_imnr * vmr_fact * dry_fact)
                else
                  scaling = scaling *          mycol_gas_imnr * vmr_fact * dry_fact
                endif
              endif
            endif

            !
            ! Interpolation of absorption coefficient and calculation of optical depth
            !
            ! Which gpoint range does this minor gas affect?
            gptS = minor_limits_gpt(1,imnr)
            gptE = minor_limits_gpt(2,imnr)

            ! Find the actual g-point to work on
            igpt = igpt0 + gptS

            ! Proceed only if the g-point is within the correct range
            if (igpt <= gptE) then
              ! What is the starting point in the stored array of minor absorption coefficients?
              minor_start = kminor_start(imnr)

              tau_minor = 0._wp
              iflav = gpt_flv(idx_tropo,igpt) ! eta interpolation depends on flavor
              minor_loc = minor_start + (igpt - gptS) ! add offset to starting point
              kminor_loc = interpolate2D(fminor(:,:,iflav,icol,ilay), kminor, minor_loc, &
                                          jeta(:,iflav,icol,ilay), myjtemp, nminork, neta, ntemp)
              tau_minor = kminor_loc * scaling

              !$acc atomic update
              tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_minor
            endif

          enddo

        enddo
      enddo
    enddo

  end subroutine gas_optical_depths_minor
  ! ----------------------------------------------------------
  !
  ! compute Rayleigh scattering optical depths
  !
  subroutine compute_tau_rayleigh(ncol,nlay,nbnd,ngpt,         &
                                  ngas,nflav,neta,npres,ntemp, &
                                  gpoint_flavor,band_lims_gpt, &
                                  krayl,                       &
                                  idx_h2o, col_dry,col_gas,    &
                                  fminor,jeta,tropo,jtemp,     &
                                  tau_rayleigh) bind(C, name="compute_tau_rayleigh")
    integer,                                     intent(in ) :: ncol,nlay,nbnd,ngpt
    integer,                                     intent(in ) :: ngas,nflav,neta,npres,ntemp
    integer,     dimension(2,ngpt),              intent(in ) :: gpoint_flavor
    integer,     dimension(2,nbnd),              intent(in ) :: band_lims_gpt ! start and end g-point for each band
    real(wp),    dimension(ngpt,neta,ntemp,2),   intent(in ) :: krayl
    integer,                                     intent(in ) :: idx_h2o
    real(wp),    dimension(ncol,nlay),           intent(in ) :: col_dry
    real(wp),    dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
    real(wp),    dimension(2,2,nflav,ncol,nlay), intent(in ) :: fminor
    integer,     dimension(2,  nflav,ncol,nlay), intent(in ) :: jeta
    logical(wl), dimension(ncol,nlay),           intent(in ) :: tropo
    integer,     dimension(ncol,nlay),           intent(in ) :: jtemp
    ! outputs
    real(wp),    dimension(ngpt,nlay,ncol),      intent(out) :: tau_rayleigh
    ! -----------------
    ! local variables
    real(wp) :: k ! rayleigh scattering coefficient
    integer  :: icol, ilay, iflav, igpt
    integer  :: itropo
    ! -----------------

    !$acc parallel loop collapse(3)
    do ilay = 1, nlay
      do icol = 1, ncol
        do igpt = 1, ngpt
          itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
          iflav = gpoint_flavor(itropo, igpt)
          k = interpolate2D(fminor(:,:,iflav,icol,ilay), &
                            krayl(:,:,:,itropo),      &
                            igpt, jeta(:,iflav,icol,ilay), jtemp(icol,ilay), ngpt, neta, ntemp)
          tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay))
        end do
      end do
    end do
  end subroutine compute_tau_rayleigh

  ! ----------------------------------------------------------
  !
  ! Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
  !   using Rayleigh scattering phase function
  !
  subroutine combine_and_reorder_nstr(ncol, nlay, ngpt, nmom, tau_abs, tau_rayleigh, tau, ssa, p) &
      bind(C, name="combine_and_reorder_nstr")
    integer, intent(in) :: ncol, nlay, ngpt, nmom
    real(wp), dimension(ngpt,nlay,ncol), intent(in ) :: tau_abs, tau_rayleigh
    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa
    real(wp), dimension(ncol,nlay,ngpt,nmom), &
                                         intent(inout) :: p
    ! -----------------------
    integer :: icol, ilay, igpt, imom
    real(wp) :: t
    ! -----------------------
    !$acc parallel loop collapse(3) &
    !$acc&     copy(tau, ssa, p) &
    !$acc&     copyin(tau_rayleigh(:ngpt,:nlay,:ncol),tau_abs(:ngpt,:nlay,:ncol))
    do icol = 1, ncol
      do ilay = 1, nlay
        do igpt = 1, ngpt
          t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
          tau(icol,ilay,igpt) = t
          if(t > 2._wp * tiny(t)) then
            ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t
          else
            ssa(icol,ilay,igpt) = 0._wp
          end if
          do imom = 1, nmom
            p(imom,icol,ilay,igpt) = 0.0_wp
          end do
          if(nmom >= 2) p(2,icol,ilay,igpt) = 0.1_wp
        end do
      end do
    end do
  end subroutine combine_and_reorder_nstr
end module mo_gas_optics_kernels
