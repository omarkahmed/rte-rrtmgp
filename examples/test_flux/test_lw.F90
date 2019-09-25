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
! Example program to demonstrate the calculation of longwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_lw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
!
! -------------------------------------------------------------------------------------------------
!
! Error checking: Procedures in rte+rrtmgp return strings which are empty if no errors occured
!   Check the incoming string, print it out and stop execution if non-empty
!
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "test_flux_compute stopping"
    stop
  end if
end subroutine stop_on_err

program test_lw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations would use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
#undef CLEAR_SKY
#undef TWO_STREAM
! #define CLEAR_SKY
! #define TWO_STREAM
  use mo_cloud_optics,            only: ty_cloud_optics
  use mo_load_cloud_coefficients,     only: load_cld_lutcoeff, load_cld_padecoeff, &
                                            read_cldpp, write_cldop, is_lw

  use mo_optical_props,       only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_optical_props_add,   only: ty_optical_props_tip
  !
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly...
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! ... and another type to encapsulate the longwave source functions.
  !
  use mo_source_functions,   only: ty_source_func_lw
  !
  ! RTE longwave driver
  !
  use mo_rte_lw,             only: rte_lw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  !
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  !
  ! modules for reading and writing files
  !
  ! RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  !
  use mo_load_coefficients,  only: load_and_init
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_lw_bc, determine_gas_names


  use mo_test_files_io,     only: read_atmos, read_direction, read_lw_bc, read_lw_rt, read_lw_Planck_sources

  implicit none
#ifdef USE_TIMING
  !
  ! Timing library
  !
#include "gptl.inc"
  ! use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, gptlfinalize, gptlsetoption, &
                                   ! gptlpercent, gptloverhead
#endif
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=132) :: flxdn_file, flxup_file
  integer            :: ncol, nlay, nbnd, n_quad_angles = 2
  logical            :: top_at_1
  integer            :: icol

  real(wp), dimension(:,:)  , target, allocatable :: flux_up, flux_dn
  real(wp), dimension(: ),            allocatable :: sfc_t     ! block_size, nblocks (emissivity is spectrally constant)
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis  ! block_size, nblocks (emissivity is spectrally constant)

  integer                                         :: nband, ngpt, nrghice

  !
  ! Classes used by rte+rrtmgp
  !
  type(ty_gas_optics_rrtmgp)  :: k_dist
  type(ty_source_func_lw)     :: sourceClr
  type(ty_source_func_lw)     :: source2str
  type(ty_source_func_lw)     :: sourceTip
  type(ty_source_func_lw)     :: source

  
  type(ty_fluxes_broadband)   :: fluxes
  type(ty_fluxes_broadband)   :: fluxesClr
  type(ty_fluxes_broadband)   :: fluxes2str
  type(ty_fluxes_broadband)   :: fluxesTip

  real(wp), dimension(:,:)  , target, allocatable :: flux_upC, flux_dnC
  real(wp), dimension(:,:)  , target, allocatable :: flux_up2, flux_dn2
  real(wp), dimension(:,:)  , target, allocatable :: flux_upP, flux_dnP
  !
  ! ty_gas_concentration holds multiple columns; we make an array of these objects to
  !   leverage what we know about the input file
  !
  class(ty_optical_props_arry), allocatable :: optical_props
  character(len=256) :: input_file

  real(wp), dimension(:,:),   allocatable :: col_dry
  type(ty_gas_concs)                      :: gas_concs
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev, t_lev

  ! clouds 
  type(ty_cloud_optics)                     :: cloud_spec
  character(len=128)                :: cloud_coeff_file = 'rrtmgp-lw-inputs-cloud-optics.nc'

  real(wp), dimension(:,:), allocatable :: cldfrac    ! cloud fraction
  real(wp), dimension(:,:), allocatable :: ciwp       ! cloud ice water path
  real(wp), dimension(:,:), allocatable :: clwp       ! cloud liquid water path
  real(wp), dimension(:,:), allocatable :: rei        ! cloud ice particle effective size (microns)
  real(wp), dimension(:,:), allocatable :: rel        ! cloud liquid particle effective radius (microns)
  integer :: icergh
  class(ty_optical_props_2str), allocatable :: cloud_optical_props
  real(wp), parameter                   :: mask_min_value = 1.e-20 ! Example; user should set as desired.
  logical,  dimension(:,:), allocatable :: cld_mask_liq, & ! Masks for liquid and ice clouds
                                           cld_mask_ice

  real(wp), dimension(:,:), allocatable :: disDFLUX
  real(wp), dimension(:,:), allocatable :: disUFLUX
  real(wp), dimension(:,:), allocatable :: tangDFLUX
  real(wp), dimension(:,:), allocatable :: tangUFLUX
  real(wp) :: dlwc
  class(ty_optical_props_1scl), allocatable :: optPropClr
  class(ty_optical_props_2str), allocatable :: optProp2str
  class(ty_optical_props_tip),  allocatable :: optPropTip
  integer                                   :: sLev, eLev
  real(wp), dimension(4)                    :: lC, iC

  logical                                   :: noSurfaceSource
  logical                                   :: isothermalAtmosphere
  logical                                   :: noLiquidCloud
  logical                                   :: noIceCloud
  integer                                   :: cloudModel
#ifdef USE_TIMING
  integer :: ret, i, jj, igas
#endif


  noSurfaceSource = .true.
  isothermalAtmosphere= .true.
  noLiquidCloud = .false.
  noIceCloud = .false.
  cloudModel = 1

  input_file="rrtmgp-lw-inputs-outputs-cloud.nc"
  input_file='/OSS_Scratch/RRTMGP/rte-rrtmgp-test-val/test/cloud_optics/ref/rrtmgp-lw-inputs-outputs-cloud.nc'
  call load_cld_lutcoeff (cloud_spec, cloud_coeff_file)
  nbnd    = cloud_spec%get_nband()
  nrghice = cloud_spec%get_num_ice_roughness_types()

  call load_cld_padecoeff(cloud_spec, cloud_coeff_file)

  print *,'nbnd', nbnd
  print *,'nrghice', nrghice

  call read_cldpp(input_file, cldfrac, clwp, ciwp, rei, rel, icergh)
  call stop_on_err(cloud_spec%set_ice_roughness(icergh))


  call read_atmos(input_file, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry)
  call load_and_init(k_dist,'coefficients_lw.nc', gas_concs)

  ncol = size(p_lay, 1)
  nlay = size(p_lay, 2)
  ngpt = k_dist%get_ngpt()
  nbnd = k_dist%get_nband()

  print *,'ncol, nlay, ngpt, nbnd', ncol, nlay, ngpt, nbnd
  print *,''
  allocate(disUFLUX(0:nlay,ncol), disDFLUX(0:nlay,ncol))
  allocate(tangUFLUX(0:nlay,ncol), tangDFLUX(0:nlay,ncol))

#ifdef CLEAR_SKY
  allocate(ty_optical_props_1scl::optical_props)
#elif defined(TWO_STREAM)
  allocate(ty_optical_props_2str::optical_props)
#else  
  allocate(ty_optical_props_tip::optical_props)
#endif
  allocate(cloud_optical_props)

if(k_dist%source_is_internal()) then
  print *, " Compute longwave gas optical depths"
  call stop_on_err(source    %alloc(ncol, nlay, k_dist))
  call stop_on_err(sourceClr %alloc(ncol, nlay, k_dist))
  call stop_on_err(source2str%alloc(ncol, nlay, k_dist))
  call stop_on_err(sourceTip %alloc(ncol, nlay, k_dist))


  call read_lw_bc     (input_file, sfc_t, sfc_emis)

else
  print *, "ERROR: Compute shortwave gas optical depths"
  STOP
end if

  !
  ! Optical properties arrays
  !
  allocate(ty_optical_props_1scl::optPropClr)
  allocate(ty_optical_props_2str::optProp2str)
  allocate(ty_optical_props_tip ::optPropTip)

  call stop_on_err(       optical_props%init(k_dist))
  call stop_on_err(       optPropClr %init(k_dist))
  call stop_on_err(       optProp2str%init(k_dist))
  call stop_on_err(       optPropTip %init(k_dist))

  call stop_on_err(optPropClr%alloc_1scl(ncol, nlay))
  call stop_on_err(optProp2str%alloc_2str(ncol, nlay))
  call stop_on_err(optPropTip %alloc_2str(ncol, nlay))

  select type (optical_props)
    type is (ty_optical_props_1scl) ! two-stream calculation
        call stop_on_err(optical_props%alloc_1scl(ncol, nlay))
    type is (ty_optical_props_2str) ! two-stream calculation
        call stop_on_err(optical_props%alloc_2str(ncol, nlay))
    type is (ty_optical_props_tip)
        call stop_on_err(optical_props%alloc_2str(ncol, nlay))
  end select

  call stop_on_err(cloud_optical_props%alloc_2str(ncol, nlay, cloud_spec))
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

  if (top_at_1) then
    sLev=1
    do sLev = 1, ncol-1
     if ( p_lev(1, sLev)/1e2 >= 6.633239E+02 ) then
        exit
      endif
    enddo
    sLev = sLev -1
    do eLev = sLev, ncol-1
     if ( p_lev(1, eLev)/1e2 >= 6.633239E+02 +  1.397207E+02 ) then
        exit
      endif
    enddo
  else
    if (cloudModel==1) then
      ! MODEL 1
      lC = [6.633239E+02,  1.397207E+02 , 3.826371E-02,  1.301899E+01]
      iC = [3.974811E+02,  9.594205E+01 , 4.415885E-02, 7.351759E+01]
    else
      ! MODEL 2
      lC = [7.439591E+02,  1.225180E+02,  4.527680E-02,  3.152708E+01]
      iC = [4.005241E+02,  9.172126E+01,  2.947712E-02,  5.827560E+01]
    endif
      
    eLev=1
    do eLev = 1, ncol-1
     if ( p_lev(1, eLev)/1e2 < lC(1) +lC(2)) then
        exit
      endif
    enddo
    eLev = eLev -1
    do sLev = eLev, ncol-1
     if ( p_lev(1, sLev)/1e2 <  lC(1)  ) then
        exit
      endif
    enddo
    clwp(1,:) = 0.
    ! rel(1,:)  = 1.
    dlwc =  lC(3)* 1e3 /(p_lev(1, eLev)-p_lev(1, sLev)) 
    print *, 'dlwc', dlwc
    do jj = eLev, sLev-1
      cldfrac(1,jj) = 1.
      clwp(1,jj) = p_lay(1,jj) * dlwc 
      rel (1,jj) = lC(4)
      print *, clwp(1,jj), p_lay(1,jj), p_lev(1,jj) -  p_lev(1,jj+1)
    enddo

    eLev=1
    do eLev = 1, ncol-1
     if ( p_lev(1, eLev)/1e2 < iC(1) +iC(2)) then
        exit
      endif
    enddo
    eLev = eLev -1
    do sLev = eLev, ncol-1
     if ( p_lev(1, sLev)/1e2 <  iC(1)  ) then
        exit
      endif
    enddo
    ciwp(1,:) = 0.
    ! rel(1,:)  = 1.
    dlwc =  iC(3)* 1e3 /(p_lev(1, eLev)-p_lev(1, sLev)) 
    print *, 'diwc', dlwc
    do jj = eLev, sLev-1
      cldfrac(1,jj) = 1.
      ciwp(1,jj) = p_lay(1,jj) * dlwc 
      rei (1,jj) = iC(4)
      print *, ciwp(1,jj), p_lay(1,jj), p_lev(1,jj) -  p_lev(1,jj+1)
    enddo


  endif

! reverse
  ! top_at_1 = .not. top_at_1
  ! do jj=1,ncol
  !   call flip_ud(nlay+1, p_lev(jj,:))
  !   call flip_ud(nlay+1, t_lev(jj,:))
  !   call flip_ud(nlay,   col_dry(jj,:))
  !   call flip_ud(nlay,   t_lay(jj,:))
  !   call flip_ud(nlay,   p_lay(jj,:))
  !   call flip_ud(nlay,   ciwp(jj,:))
  !   call flip_ud(nlay,   rei(jj,:))
  !   call flip_ud(nlay,   clwp(jj,:))
  !   call flip_ud(nlay,   rel(jj,:))
  ! enddo
  ! gas

  do igas=1,size(gas_concs%concs)
  do jj=1,ncol
    call flip_ud(nlay, gas_concs%concs(igas)%conc(jj,:))
  enddo
  enddo


  n_quad_angles=4
  ! possibility to zero out liquid or water cloud component
  if (noLiquidCloud) clwp = 0.
  if (noIceCloud)    ciwp=0.

  ! set isothermal atmosphere
  if (isothermalAtmosphere) then
    t_lev=300.
    t_lay=300.
    sfc_t=300.
  endif
    
  cld_mask_liq = (cldfrac >= mask_min_value .and. clwp >= mask_min_value)
  cld_mask_ice = (cldfrac >= mask_min_value .and. ciwp >= mask_min_value)


  call stop_on_err(cloud_spec%cloud_optics(ncol, nlay, nbnd, nrghice,  &
                                           cld_mask_liq, cld_mask_ice, &
                                           clwp,         ciwp,         &
                                           rel,          rei, &
                                           cloud_optical_props) )  

  ! Number of quadrature angles
  ! call read_lw_rt(input_file, n_quad_angles)
  print *,'---------------------- cloud'
  print *,'nlay', nlay
  print *,'n_quad_angles', n_quad_angles
  print *,'TOP:', top_at_1, p_lay(1, 1), p_lay(1, nlay)

  call stop_on_err(k_dist%gas_optics_int_cloud(p_lay, p_lev, &
                                        t_lay, sfc_t, &
                                        gas_concs,    &
                                        optical_props,&
                                        source,     &
                                        tlev=t_lev, &
                                        col_dry=col_dry,&
                                        optical_props_clouds=cloud_optical_props) )

  call stop_on_err(k_dist%gas_optics_int_cloud(p_lay, p_lev, &
                                        t_lay, sfc_t, &
                                        gas_concs,  &
                                        optPropClr, &
                                        sourceClr,  &
                                        tlev=t_lev, &
                                        col_dry=col_dry,&
                                        optical_props_clouds=cloud_optical_props) )
  call stop_on_err(k_dist%gas_optics_int_cloud(p_lay, p_lev, &
                                        t_lay, sfc_t, &
                                        gas_concs,                       &
                                        optProp2str,                   &
                                        source2str,                         &
                                        tlev=t_lev, &
                                        col_dry=col_dry,&
                                        optical_props_clouds=cloud_optical_props) )
  call stop_on_err(k_dist%gas_optics_int_cloud(p_lay, p_lev, &
                                        t_lay, sfc_t, &
                                        gas_concs,                       &
                                        optPropTip,                   &
                                        sourceTip,                         &
                                        tlev=t_lev, &
                                        col_dry=col_dry,&
                                        optical_props_clouds=cloud_optical_props) )
  call stop_on_err(optPropTip%delta_scale())


! zero out the surface source
if (noSurfaceSource) then
  source    %sfc_source=0.
  sourceClr %sfc_source=0.
  source2str%sfc_source=0.
  sourceTip %sfc_source=0.
  source2str%planckSfc=0.
endif

  allocate(    flux_up (ncol,nlay+1     ),     flux_dn (ncol,nlay+1     ))
  allocate(    flux_upC(ncol,nlay+1     ),     flux_dnC(ncol,nlay+1     ))
  allocate(    flux_up2(ncol,nlay+1     ),     flux_dn2(ncol,nlay+1     ))
  allocate(    flux_upP(ncol,nlay+1     ),     flux_dnP(ncol,nlay+1     ))

  fluxes%flux_up     => flux_up
  fluxes%flux_dn     => flux_dn
  fluxesClr%flux_up  => flux_upC
  fluxesClr%flux_dn  => flux_dnC
  fluxes2str%flux_up => flux_up2
  fluxes2str%flux_dn => flux_dn2
  fluxesTip %flux_up => flux_upP
  fluxesTip %flux_dn => flux_dnP


#ifdef USE_TIMING
  !
  ! Initialize timers
  !
  ret = gptlsetoption (gptlpercent, 1)        ! Turn on "% of" print
  ret = gptlsetoption (gptloverhead, 0)       ! Turn off overhead estimate
  ret = gptlinitialize()
#endif
  print *,'Optical properties'


#ifdef USE_TIMING
  do i = 1, 10

      ret =  gptlstart('rte_lw_clr')
#endif

  call stop_on_err(rte_lw(optPropClr,             &
                          top_at_1,              &
                          sourceClr,        &
                          sfc_emis, &
                          fluxesClr, n_gauss_angles = n_quad_angles))

#ifdef USE_TIMING
    ret =  gptlstop('rte_lw_clr')

    ret =  gptlstart('rte_lw_2str')
#endif

  call stop_on_err(rte_lw(optProp2str,             &
                          top_at_1,              &
                          source2str,        &
                          sfc_emis, &
                          fluxes2str))

#ifdef USE_TIMING
    ret =  gptlstop('rte_lw_2str')

    ret =  gptlstart('rte_lw_tip')
#endif

  call stop_on_err(rte_lw(optPropTip,             &
                          top_at_1,              &
                          sourceTip,        &
                          sfc_emis, &
                          fluxesTip, n_gauss_angles = n_quad_angles))

#ifdef USE_TIMING
    ret =  gptlstop('rte_lw_tip')

    ret =  gptlstart('rte_lw_ttang')
#endif

  call RTREGCLDIP( ncol, nlay, ngpt, n_quad_angles, &
               optPropTip, &
               sfc_emis,&
               source2str%planckFrac,&
               source2str%planckLay,&
               source2str%planckLev,&
               source2str%planckSfc,&
               source2str%lay_source, &
               source2str%lev_source_inc, &
               source2str%lev_source_dec,&
               tangUFLUX, tangDFLUX, top_at_1)
  ! End timers
  !
#ifdef USE_TIMING
    ret =  gptlstop('rte_lw_ttang')
  end do

    ret =  gptlstart('rte_disort')
#endif

! zero out the surface source
  if (noSurfaceSource) then
     sfc_t=1e-6
  endif
  call RTRDIS( ncol, nlay, ngpt, optProp2str, &
               t_lev, sfc_t,sfc_emis,&
               source2str%planckFrac,&
               source2str%planckLev,&
               disUFLUX, disDFLUX, top_at_1)
  ! End timers
  !
#ifdef USE_TIMING
    ret =  gptlstop('rte_disort')

  !
  ! End timers
  !
  ret=gptlpr_file('timing')
  ret = gptlfinalize()
#endif

do icol=1,1!ncol  
print *,'UP FLUX'
print *,'  clear     2 str     Tang IP   DISORT      clear     2str      tangip'
do jj=1,nlay+1
  if (top_at_1) then
    i = nlay - jj +1
  else
    i = jj-1
  endif
  print '(I4, 100F10.4)', jj, flux_upC(icol,jj), flux_up2(icol,jj), &
     flux_upP(icol,jj),  disUFLUX(i,icol),&
     flux_upC(icol,jj)-disUFLUX(i,icol), flux_up2(icol,jj)-disUFLUX(i,icol),&
     flux_upP(icol,jj)-disUFLUX(i,icol),&
     tangUFLUX(i,icol), tangUFLUX(i,icol)-disUFLUX(i,icol)
enddo  
print *,'------- '
print *,'DOWN FLUX'
print *,'  clear     2 str     Tang IP   DISORT      clear     2str      tangip'
 
do jj=1,nlay+1
  if (top_at_1) then
    i = nlay - jj +1
  else
    i = jj-1
  endif
  print '(I4, 100F10.4)', jj, flux_dnC(icol,jj), flux_dn2(icol,jj), &
     flux_dnP(icol,jj), disDFLUX(i,icol),&
     flux_dnC(icol,jj)- disDFLUX(i,icol), flux_dn2(icol,jj)- disDFLUX(i,icol),&
     flux_dnP(icol,jj)- disDFLUX(i,icol),&
     tangDFLUX(i,icol), tangDFLUX(i,icol)-disDFLUX(i,icol)
enddo  
print *,'------- '
do jj=1,nlay+1
  if (top_at_1) then
    i = nlay - jj +1
  else
    i = jj-1
  endif
  write(1001, *) flux_upC(icol,jj), flux_up2(icol,jj), &
     flux_upP(icol,jj), disUFLUX(i,icol),&
     flux_upC(icol,jj)-disUFLUX(i,icol), flux_up2(icol,jj)-disUFLUX(i,icol),&
     flux_upP(icol,jj)-disUFLUX(i,icol),&
     tangUFLUX(i,icol), tangUFLUX(i,icol)-disUFLUX(i,icol)
enddo  
do jj=1,nlay+1
  if (top_at_1) then
    i = nlay - jj +1
  else
    i = jj-1
  endif
  write(1002, *) flux_dnC(icol,jj), flux_dn2(icol,jj), &
     flux_dnP(icol,jj), disDFLUX(i,icol),&
     flux_dnC(icol,jj)- disDFLUX(i,icol), flux_dn2(icol,jj)- disDFLUX(i,icol),&
     flux_dnP(icol,jj)- disDFLUX(i,icol),&
     tangDFLUX(i,icol), tangDFLUX(i,icol)-disDFLUX(i,icol)
  write(1003, *) p_lev(1,jj)
enddo  
enddo
do jj=1,nlay+1
   write(1003, *) p_lev(1,jj)
enddo  
print *,'cloud OT', sum(cloud_optical_props%tau(1,:,1))
contains
  subroutine flip_ud(nlay, dat)
    integer,                   intent(in)    :: nlay
    real(wp), dimension(nlay), intent(inout) :: dat
    real(wp), dimension(nlay) :: buf
    buf=dat
    dat=buf(nlay:1:-1)
  end subroutine flip_ud
    
end program test_lw














