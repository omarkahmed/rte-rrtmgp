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
! Kernels for computing broadband fluxes by summing over all elements in the spectral dimension
!
! -------------------------------------------------------------------------------------------------
module mo_fluxes_broadband_kernels
  use, intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp
  implicit none
  private
  public :: sum_broadband, net_broadband


  interface

    subroutine sum_broadband(ncol, nlev, ngpt, spectral_flux, broadband_flux) bind(C, name="sum_broadband")
      use mo_rte_kind, only: wp
      integer, value,                        intent(in ) :: ncol, nlev, ngpt
      real(wp), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux
      real(wp), dimension(ncol, nlev),       intent(out) :: broadband_flux
    end subroutine

    subroutine net_broadband_full(ncol, nlev, ngpt, spectral_flux_dn, spectral_flux_up, broadband_flux_net) &
      bind(C, name="net_broadband_full")
      use mo_rte_kind, only: wp
      integer, value,                        intent(in ) :: ncol, nlev, ngpt
      real(wp), dimension(ncol, nlev, ngpt), intent(in ) :: spectral_flux_dn, spectral_flux_up
      real(wp), dimension(ncol, nlev),       intent(out) :: broadband_flux_net
    end subroutine

    subroutine net_broadband_precalc(ncol, nlev, flux_dn, flux_up, broadband_flux_net) &
      bind(C, name="net_broadband_precalc")
      use mo_rte_kind, only: wp
      integer, value,                  intent(in ) :: ncol, nlev
      real(wp), dimension(ncol, nlev), intent(in ) :: flux_dn, flux_up
      real(wp), dimension(ncol, nlev), intent(out) :: broadband_flux_net
    end subroutine

  end interface


  interface net_broadband
    procedure net_broadband_precalc
    procedure net_broadband_full
  end interface net_broadband


end module mo_fluxes_broadband_kernels
