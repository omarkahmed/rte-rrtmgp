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
! Encapsulate optical properties defined on a spectral grid of N bands.
!   The bands are described by their limiting wavenumbers. They need not be contiguous or complete.
!   A band may contain more than one spectral sub-point (g-point) in which case a mapping must be supplied.
!   A name may be provided and will be prepended to error messages.
!   The base class (ty_optical_props) encapsulates only this spectral discretization and must be initialized
!      with the spectral information before use.
!
!   Optical properties may be represented as arrays with dimensions ncol, nlay, ngpt
!   (abstract class ty_optical_props_arry).
!   The type holds arrays depending on how much information is needed
!   There are three possibilites
!      ty_optical_props_1scl holds absorption optical depth tau, used in calculations accounting for extinction and emission
!      ty_optical_props_2str holds extincion optical depth tau, single-scattering albedo ssa, and
!        asymmetry parameter g. These fields are what's needed for two-stream calculations.
!      ty_optical_props_nstr holds extincion optical depth tau, single-scattering albedo ssa, and
!        phase function moments p with leading dimension nmom. These fields are what's needed for multi-stream calculations.
!   These classes must be allocated before use. Initialization and allocation can be combined.
!   The classes have a validate() function that checks all arrays for valid values (e.g. tau > 0.)
!
! Optical properties can be delta-scaled (though this is currently implemented only for two-stream arrays)
!
! Optical properties can increment or "add themselves to" a set of properties represented with arrays
!   as long as both sets have the same underlying band structure. Properties defined by band
!   may be added to properties defined by g-point; the same value is assumed for all g-points with each band.
!
! Subsets of optical properties held as arrays may be extracted along the column dimension.
!
! -------------------------------------------------------------------------------------------------
module mo_optical_props_add
  use mo_rte_kind,              only: wp
  use mo_optical_props,         only: ty_optical_props_2str
 
  ! Implemented based on the paper
  ! Tang G, P Yang, GW Kattawar, X Huang, EJ Mlawer, BA Baum, MD King, 2018: Improvement of 
  ! the Simulation of Cloud Longwave Scattering in Broadband Radiative Transfer Models, 
  ! Journal of the Atmospheric Sciences 75 (7), 2217-2233
  ! https://doi.org/10.1175/JAS-D-18-0014.1
  !
  type, extends(ty_optical_props_2str) :: ty_optical_props_tip
  contains
      procedure, public  :: delta_scale => delta_scale_tip
  end type  ty_optical_props_tip


contains  
  ! ------------------------------------------------------------------------------------------
  ! --- delta scaling
  ! ------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------
  function delta_scale_tip(this, for) result(err_message)
    class(ty_optical_props_tip), intent(inout) :: this
    real(wp), dimension(:,:,:), optional, &
                                  intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    real(wp) :: wf
    integer  :: icol, ilay, igpt
    integer :: ncol, nlay, ngpt
    ! --------------------------------
    ncol = this%get_ncol()
    nlay = this%get_nlay()
    ngpt = this%get_ngpt()
    err_message = ""
    call scalingTang(ncol, nlay, ngpt, this%tau, this%ssa, this%g)
  end function delta_scale_tip  
  

  pure subroutine scalingTang(ncol, nlay, ngpt, tau, ssa, g)
    integer ,                              intent(in)    :: ncol
    integer ,                              intent(in)    :: nlay
    integer ,                              intent(in)    :: ngpt
    real(wp), dimension(ncol, nlay, ngpt), intent(inout) :: tau
    real(wp), dimension(ncol, nlay, ngpt), intent(inout) :: ssa
    real(wp), dimension(ncol, nlay, ngpt), intent(inout) :: g

    integer  :: icol, ilay, igpt
    real(wp) :: gl, ssal
    do igpt=1,ngpt
      do ilay=1,nlay
        do icol=1,ncol
          gl = (1._wp + g(icol, ilay, igpt)) / 2._wp
          ssal = ssa(icol, ilay, igpt)
          ! Eq.15 of the paper
          tau(icol, ilay, igpt) = (1._wp - ssal * gl) * tau(icol, ilay, igpt)
          ! 
          ! here ssa is used to store parameter wb/[(]1-w(1-b)] of Eq.21 of the Tang's paper
          ! actually it is in line of parameter rescaling defined in Eq.7
          !
          ! here it is a good place to add factor 0.5 (0.4 note A of Table) to save a flop
          !
          ssa(icol, ilay, igpt) = (1._wp - gl) *ssal / (1._wp - ssal * gl)*0.4/0.5
          g(icol, ilay, igpt)   = gl
        enddo
      enddo
    enddo
  end subroutine scalingTang
end module mo_optical_props_add
