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
!
! The optional solar variability extension class used by RRMTGP must be initialized with data stored 
! in a netCDF file.
! RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
! local environment. This module provides a straight-forward implementation of reading the data.
!
! -------------------------------------------------------------------------------------------------
module mo_load_solar_var_tables
  !
  ! Modules for working with rte and rrtmgp
  !
  use mo_rte_kind,           only: wp, wl
  use mo_solar_variability,  only: ty_solar_var

  ! --------------------------------------------------
  use mo_simple_netcdf,      only: read_field, get_dim_size
  use netcdf

  implicit none
  private
  public :: load_solar_var_tables

contains
  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
  ! read solar variabilty tables from NetCDF file
  subroutine load_solar_var_tables(solar_var, filename)
    class(ty_solar_var),  intent(inout) :: solar_var ! solar variability data
    character(len=*),     intent(in   ) :: filename

    real(wp), dimension(:,:    ), allocatable :: svar_avgcyc
    ! -----------------
    !
    ! Book-keeping variables
    !
    integer :: ncid
    integer :: nsolarterms,     &
               nsolarfrac
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_solar_var_tables(): can't open file " // trim(fileName))
    nsolarterms         = get_dim_size(ncid,'n_solar_terms')
    nsolarfrac          = get_dim_size(ncid,'n_solar_frac')
    ! -----------------
    !
    ! Read the arrays
    !
    allocate (svar_avgcyc(nsolarterms,nsolarfrac))

    svar_avgcyc = read_field(ncid, 'solar_var_avgcyc', nsolarterms, nsolarfrac)

    call stop_on_err(solar_var%load_avgcyc(svar_avgcyc))

    if (allocated(svar_avgcyc)) deallocate(svar_avgcyc)

    ! --------------------------------------------------
    ncid = nf90_close(ncid)
  end subroutine load_solar_var_tables

end module mo_load_solar_var_tables
