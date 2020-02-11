! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Kernels to permute arrays

module mo_rrtmgp_util_reorder_kernels
  use mo_rte_kind,      only: wp
  implicit none

  interface
    subroutine reorder_123x321_kernel(d1, d2, d3, array_in, array_out) bind(C, name="reorder_123x321_kernel")
      use mo_rte_kind,      only: wp
      implicit none
      integer, value                 , intent( in) :: d1, d2, d3
      real(wp), dimension(d1, d2, d3), intent( in) :: array_in
      real(wp), dimension(d3, d2, d1), intent(out) :: array_out
    end subroutine

    subroutine reorder_123x312_kernel(d1, d2, d3, array_in, array_out) bind(C, name = "reorder_123x312_kernel")
      use mo_rte_kind,      only: wp
      integer, value,                  intent( in) :: d1, d2, d3
      real(wp), dimension(d1, d2, d3), intent( in) :: array_in
      real(wp), dimension(d3, d1, d2), intent(out) :: array_out
    end subroutine
  end interface

end module mo_rrtmgp_util_reorder_kernels
