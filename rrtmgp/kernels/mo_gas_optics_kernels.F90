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

    subroutine compute_tau_rayleigh(ncol,nlay,nbnd,ngpt,         &
                                    ngas,nflav,neta,npres,ntemp, &
                                    gpoint_flavor,band_lims_gpt, &
                                    krayl,                       &
                                    idx_h2o, col_dry,col_gas,    &
                                    fminor,jeta,tropo,jtemp,     &
                                    tau_rayleigh) bind(C, name="compute_tau_rayleigh")
      use mo_rte_kind,      only : wp, wl
      implicit none
      integer, value,                              intent(in ) :: ncol,nlay,nbnd,ngpt
      integer, value,                              intent(in ) :: ngas,nflav,neta,npres,ntemp
      integer,     dimension(2,ngpt),              intent(in ) :: gpoint_flavor
      integer,     dimension(2,nbnd),              intent(in ) :: band_lims_gpt ! start and end g-point for each band
      real(wp),    dimension(ngpt,neta,ntemp,2),   intent(in ) :: krayl
      integer, value,                              intent(in ) :: idx_h2o
      real(wp),    dimension(ncol,nlay),           intent(in ) :: col_dry
      real(wp),    dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
      real(wp),    dimension(2,2,nflav,ncol,nlay), intent(in ) :: fminor
      integer,     dimension(2,  nflav,ncol,nlay), intent(in ) :: jeta
      logical(wl), dimension(ncol,nlay),           intent(in ) :: tropo
      integer,     dimension(ncol,nlay),           intent(in ) :: jtemp
      real(wp),    dimension(ngpt,nlay,ncol),      intent(out) :: tau_rayleigh
    end subroutine

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
      use mo_rte_kind,      only : wp, wl
      implicit none
      integer, value,                        intent(in) :: ncol,nlay,nbnd,ngpt
      integer, value,                        intent(in) :: ngas,nflav,neta,npres,ntemp
      integer, value,                        intent(in) :: nminorlower, nminorklower,nminorupper, nminorkupper
      integer, value,                        intent(in) :: idx_h2o
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
      real(wp), dimension(2,    nflav,ncol,nlay       ), intent(in) :: col_mix
      real(wp), dimension(2,2,2,nflav,ncol,nlay       ), intent(in) :: fmajor
      real(wp), dimension(2,2,  nflav,ncol,nlay       ), intent(in) :: fminor
      real(wp), dimension(            ncol,nlay       ), intent(in) :: play, tlay      ! pressure and temperature
      real(wp), dimension(            ncol,nlay,0:ngas), intent(in) :: col_gas
      integer,  dimension(2,    nflav,ncol,nlay       ), intent(in) :: jeta
      integer,  dimension(            ncol,nlay       ), intent(in) :: jtemp
      integer,  dimension(            ncol,nlay       ), intent(in) :: jpress
      real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    end subroutine

  end interface


contains


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
