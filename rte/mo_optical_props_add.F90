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
  ! ! -------------------------------------------------------------------------------------------------
  ! !
  ! ! Base class for optical properties
  ! !   Describes the spectral discretization including the wavenumber limits
  ! !   of each band (spectral region) and the mapping between g-points and bands
  ! !
  ! ! -------------------------------------------------------------------------------------------------
  ! type, public :: ty_optical_props
  !   integer,  dimension(:,:), allocatable :: band2gpt       ! (begin g-point, end g-point) = band2gpt(2,band)
  !   integer,  dimension(:),   allocatable :: gpt2band       ! band = gpt2band(g-point)
  !   real(wp), dimension(:,:), allocatable :: band_lims_wvn  ! (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  !   character(len=name_len)               :: name = ""
  ! contains
  !   generic,   public  :: init => init_base, init_base_from_copy
  !   procedure, private :: init_base
  !   procedure, private :: init_base_from_copy
  !   procedure, public  :: is_initialized => is_initialized_base
  !   procedure, private :: is_initialized_base
  !   procedure, public  :: finalize => finalize_base
  !   procedure, private :: finalize_base
  !   procedure, public  :: get_nband
  !   procedure, public  :: get_ngpt
  !   procedure, public  :: get_gpoint_bands
  !   procedure, public  :: convert_band2gpt
  !   procedure, public  :: convert_gpt2band
  !   procedure, public  :: get_band_lims_gpoint
  !   procedure, public  :: get_band_lims_wavenumber
  !   procedure, public  :: get_band_lims_wavelength
  !   procedure, public  :: bands_are_equal
  !   procedure, public  :: gpoints_are_equal
  !   procedure, public  :: expand
  !   procedure, public  :: set_name
  !   procedure, public  :: get_name
  ! end type
  !----------------------------------------------------------------------------------------
  !
  ! Optical properties as arrays, normally dimensioned ncol, nlay, ngpt/nbnd
  !   The abstract base class for arrays defines what procedures will be available
  !   The optical depth field is also part of the abstract base class, since
  !    any representation of values as arrays needs an optical depth field
  !
  ! -------------------------------------------------------------------------------------------------
  ! type, extends(ty_optical_props), abstract, public :: ty_optical_props_arry
  !   real(wp), dimension(:,:,:), allocatable :: tau ! optical depth (ncol, nlay, ngpt)
  ! contains
  !   procedure, public  :: get_ncol
  !   procedure, public  :: get_nlay
  !   !
  !   ! Increment another set of values
  !   !
  !   procedure, public  :: increment

  !   !
  !   ! Deferred procedures -- each must be implemented in each child class with
  !   !   arguments following the abstract interface (defined below)
  !   !
  !   procedure(validate_abstract),     deferred, public  :: validate
  !   procedure(delta_scale_abstract),  deferred, public  :: delta_scale
  !   procedure(subset_range_abstract), deferred, public  :: get_subset
  ! end type
  ! !
  ! ! Interfaces for the methods to be implemented
  ! !
  ! abstract interface
  !   !
  !   ! Validation function looks only at internal data
  !   !
  !   function validate_abstract(this) result(err_message)
  !     import ty_optical_props_arry
  !     class(ty_optical_props_arry),  intent(in) :: this
  !     character(len=128)  :: err_message
  !   end function validate_abstract

  !   !
  !   ! Delta-scaling
  !   !
  !   function delta_scale_abstract(this, for) result(err_message)
  !     import ty_optical_props_arry
  !     import wp
  !     class(ty_optical_props_arry),  intent(inout) :: this
  !     real(wp), dimension(:,:,:), optional, &
  !                                    intent(in   ) :: for
  !     ! Forward scattering fraction; g**2 if not provided
  !     character(len=128)  :: err_message
  !   end function delta_scale_abstract

  !   !
  !   ! Subsetting -- currently there are only routines with start col and count
  !   !
  !   function subset_range_abstract(full, start, n, subset) result(err_message)
  !     import ty_optical_props_arry
  !     class(ty_optical_props_arry), intent(inout) :: full
  !     integer,                      intent(in   ) :: start, n
  !     class(ty_optical_props_arry), intent(inout) :: subset
  !     character(128)                              :: err_message
  !   end function subset_range_abstract
  ! end interface
  !----------------------------------------------------------------------------------------
  !
  !   ty_optical_props_arry  includes only (extinction) optical depth
  !   Class two-stream adds arrays for single scattering albedo ssa and
  !     asymmetry parameter needed in two-stream methods
  !   Class n-stream adds arrays for single scattering albedo ssa and
  !     phase function moments (index 1 = g) for use with discrete ordinate methods
  !
  ! -------------------------------------------------------------------------------------------------

  ! Implemented based on the paper
  ! Tang G, P Yang, GW Kattawar, X Huang, EJ Mlawer, BA Baum, MD King, 2018: Improvement of 
  ! the Simulation of Cloud Longwave Scattering in Broadband Radiative Transfer Models, 
  ! Journal of the Atmospheric Sciences 75 (7), 2217-2233

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
          tau(icol, ilay, igpt) = (1._wp - ssal * gl) * tau(icol, ilay, igpt)
          !  parameter wb/[(]1-w(1-b)] to put in Eq.21 of Tang's paper
          !  here it is a good place to add factor 0.5 (0.4 note A of Table) to save a flop
          ssa(icol, ilay, igpt) = (1._wp - gl) *ssal / (1._wp - ssal * gl)*0.4/0.5
          g(icol, ilay, igpt)   = gl
        enddo
      enddo
    enddo
  end subroutine scalingTang
end module mo_optical_props_add
