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
! Routines for handling strings:
!   convert to lower case
!   does a string exist within an array of strings?
!   what is the location of a string with an array?
!
! -------------------------------------------------------------------------------------------------
module mo_rrtmgp_util_string
  implicit none

  ! List of character for case conversion
  character(len=26), parameter :: LOWER_CASE_CHARS = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: UPPER_CASE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  interface

    subroutine lower_case( input_string , output_string ) bind(C,name="lower_case")
      use iso_c_binding
      character(kind=c_char), intent(in   ) :: input_string (*)
      character(kind=c_char), intent(  out) :: output_string(*)
    end subroutine

  end interface

contains
  ! --------------------------------------------------------------------------------------
  !
  ! Is string somewhere in array?
  !
  function string_in_array(string, array)
    character(len=*),               intent(in) :: string
    character(len=*), dimension(:), intent(in) :: array
    logical                                    :: string_in_array

    integer :: i
    character(len=128) :: lc_string
    character(len=128) :: tmpstr

    string_in_array = .false.
    call char_f2c( string , lc_string )
    call lower_case( lc_string , lc_string )
    call char_c2f( lc_string , lc_string )
    do i = 1, size(array)
      call char_f2c( array(i) , tmpstr )
      call lower_case( tmpstr , tmpstr )
      call char_c2f( tmpstr , tmpstr );
      if ( trim(lc_string) == trim(tmpstr)) then
        string_in_array = .true.
        exit
      end if
    end do
  end function string_in_array
  ! --------------------------------------------------------------------------------------
  !
  ! Is string somewhere in array?
  !
  function string_loc_in_array(string, array)
    character(len=*),               intent(in) :: string
    character(len=*), dimension(:), intent(in) :: array
    integer                                    :: string_loc_in_array

    integer :: i
    character(len=128) :: lc_string
    character(len=128) :: tmpstr

    string_loc_in_array = -1
    call char_f2c( string , lc_string )
    call lower_case( lc_string , lc_string )
    call char_c2f( lc_string , lc_string )
    do i = 1, size(array)
      call char_f2c( array(i) , tmpstr )
      call lower_case( tmpstr , tmpstr )
      call char_c2f( tmpstr , tmpstr )
      if ( trim(lc_string) == trim(tmpstr) ) then
        string_loc_in_array = i
        exit
      end if
    end do
  end function string_loc_in_array
  ! --------------------------------------------------------------------------------------
  subroutine char_f2c( for , cpp )
    use iso_c_binding
    character(len=*), intent(in   ) :: for
    character(len=*), intent(  out) :: cpp
    integer flen
    flen = len_trim(for)
    cpp(1:flen+1) = for(1:flen)//C_NULL_CHAR
    cpp(flen+2:len(cpp)) = ' '
  end subroutine
  subroutine char_c2f( cpp , for )
    use iso_c_binding
    character(len=*), intent(in   ) :: cpp
    character(len=*), intent(  out) :: for
    integer :: loc
    loc = index(cpp,C_NULL_CHAR)
    for(1:loc-1) = cpp(1:loc-1)
    for(loc:len(for)) = ' '
  end subroutine
end module
