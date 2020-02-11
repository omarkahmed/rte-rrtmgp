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
! Routines for permuting arrays: here one (x,y,z) -> (z,x,y) and (x,y,z) -> (z,y,x)
!   Routines are just the front end to kernels
!
! -------------------------------------------------------------------------------------------------
module mo_rrtmgp_util_reorder
  implicit none
  interface

    subroutine reorder123x312(d1, d2, d3, array, array_out) bind(C,name="reorder123x312")
      use mo_rte_kind, only: wp
      integer, value, intent(in) :: d1, d2, d3
      real(wp), dimension(d1,d2,d3), intent(in ) :: array
      real(wp), dimension(d1,d2,d3), intent(out) :: array_out
    end subroutine

    subroutine reorder123x321(d1, d2, d3, array, array_out) bind(C,name="reorder123x321")
      use mo_rte_kind, only: wp
      integer, value, intent(in) :: d1, d2, d3
      real(wp), dimension(d1,d2,d3), intent(in ) :: array
      real(wp), dimension(d1,d2,d3), intent(out) :: array_out
    end subroutine
    
  end interface
end module
