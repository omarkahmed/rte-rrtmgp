module mo_load_cloud_coefficients
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_cloud_optics,  only: ty_cloud_optics
  use mo_simple_netcdf, only: read_field, read_string, var_exists, get_dim_size, &
                              write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_cld_lutcoeff, load_cld_padecoeff, read_cldpp, write_cldop
  public :: read_cldop_nbnd, write_cldop_ngpt, is_lw, is_sw
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read cloud optical property LUT coefficients from NetCDF file
  !
  subroutine load_cld_lutcoeff(cloud_spec, cld_coeff_file)
    class(ty_cloud_optics),     intent(inout) :: cloud_spec
    character(len=*),           intent(in   ) :: cld_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrghice, nsize_liq, nsize_ice

    real(wp), dimension(:,:), allocatable                :: band_lims_wvn
    ! Lookup table interpolation constants
    real(wp) :: radliq_lwr          ! liquid particle size lower bound for interpolation
    real(wp) :: radliq_upr          ! liquid particle size upper bound for interpolation
    real(wp) :: radliq_fac          ! constant for calculating interpolation indices for liquid
    real(wp) :: radice_lwr          ! ice particle size lower bound for interpolation
    real(wp) :: radice_upr          ! ice particle size upper bound for interpolation
    real(wp) :: radice_fac          ! constant for calculating interpolation indices for ice
    ! LUT coefficients
    real(wp), dimension(:,:),   allocatable :: lut_extliq   ! extinction: liquid
    real(wp), dimension(:,:),   allocatable :: lut_ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:),   allocatable :: lut_asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:), allocatable :: lut_extice   ! extinction: ice
    real(wp), dimension(:,:,:), allocatable :: lut_ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:), allocatable :: lut_asyice   ! asymmetry parameter: ice
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_lutcoeff(): can't open file " // trim(cld_coeff_file))

    ! Read LUT coefficient dimensions
    nband     = get_dim_size(ncid,'nband')
    nrghice   = get_dim_size(ncid,'nrghice')
    nsize_liq = get_dim_size(ncid,'nsize_liq')
    nsize_ice = get_dim_size(ncid,'nsize_ice')

    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber')

    ! Read LUT constants
    radliq_lwr = read_field(ncid, 'radliq_lwr')
    radliq_upr = read_field(ncid, 'radliq_upr')
    radliq_fac = read_field(ncid, 'radliq_fac')
    radice_lwr = read_field(ncid, 'radice_lwr')
    radice_upr = read_field(ncid, 'radice_upr')
    radice_fac = read_field(ncid, 'radice_fac')

    ! Allocate cloud property lookup table input arrays
    allocate(lut_extliq(nsize_liq, nband), &
             lut_ssaliq(nsize_liq, nband), &
             lut_asyliq(nsize_liq, nband), &
             lut_extice(nsize_ice, nband, nrghice), &
             lut_ssaice(nsize_ice, nband, nrghice), &
             lut_asyice(nsize_ice, nband, nrghice))

    ! Read LUT coefficients
    lut_extliq = read_field(ncid, 'lut_extliq',  nsize_liq, nband)
    lut_ssaliq = read_field(ncid, 'lut_ssaliq',  nsize_liq, nband)
    lut_asyliq = read_field(ncid, 'lut_asyliq',  nsize_liq, nband)
    lut_extice = read_field(ncid, 'lut_extice',  nsize_ice, nband, nrghice)
    lut_ssaice = read_field(ncid, 'lut_ssaice',  nsize_ice, nband, nrghice)
    lut_asyice = read_field(ncid, 'lut_asyice',  nsize_ice, nband, nrghice)

    ncid = nf90_close(ncid)

    call stop_on_err(cloud_spec%load(band_lims_wvn,                      &
                                     radliq_lwr, radliq_upr, radliq_fac, &
                                     radice_lwr, radice_upr, radice_fac, &
                                     lut_extliq, lut_ssaliq, lut_asyliq, &
                                     lut_extice, lut_ssaice, lut_asyice))
  end subroutine load_cld_lutcoeff
  !--------------------------------------------------------------------------------------------------------------------
  ! read cloud optical property Pade coefficients from NetCDF file
  !
  subroutine load_cld_padecoeff(cloud_spec, cld_coeff_file)
    class(ty_cloud_optics),       intent(inout) :: cloud_spec
    character(len=*),             intent(in   ) :: cld_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

    ! Spectral discretization
    real(wp), dimension(:,:), allocatable :: band_lims_wvn

    ! Pade coefficients
    real(wp), dimension(:,:,:),   allocatable :: pade_extliq   ! extinction: liquid
    real(wp), dimension(:,:,:),   allocatable :: pade_ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:,:),   allocatable :: pade_asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice   ! extinction: ice
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:,:), allocatable :: pade_asyice   ! asymmetry parameter: ice

    ! Pade particle size regime boundaries
    real(wp),  dimension(:),       allocatable :: pade_sizreg_extliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_asyliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_extice
    real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaice
    real(wp),  dimension(:),       allocatable :: pade_sizreg_asyice
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_padecoeff(): can't open file " // trim(cld_coeff_file))

    ! Read Pade coefficient dimensions
    nband        = get_dim_size(ncid,'nband')
    nrghice      = get_dim_size(ncid,'nrghice')
    nsizereg     = get_dim_size(ncid,'nsizereg')
    ncoeff_ext   = get_dim_size(ncid,'ncoeff_ext')
    ncoeff_ssa_g = get_dim_size(ncid,'ncoeff_ssa_g')
    nbound       = get_dim_size(ncid,'nbound')

    !
    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber')

    ! Allocate cloud property Pade coefficient input arrays
    allocate(pade_extliq(nband, nsizereg, ncoeff_ext),   &
             pade_ssaliq(nband, nsizereg, ncoeff_ssa_g), &
             pade_asyliq(nband, nsizereg, ncoeff_ssa_g), &
             pade_extice(nband, nsizereg, ncoeff_ext,   nrghice), &
             pade_ssaice(nband, nsizereg, ncoeff_ssa_g, nrghice), &
             pade_asyice(nband, nsizereg, ncoeff_ssa_g, nrghice))

    pade_extliq  = read_field(ncid, 'pade_extliq', nband, nsizereg, ncoeff_ext)
    pade_ssaliq  = read_field(ncid, 'pade_ssaliq', nband, nsizereg, ncoeff_ssa_g)
    pade_asyliq  = read_field(ncid, 'pade_asyliq', nband, nsizereg, ncoeff_ssa_g)
    pade_extice  = read_field(ncid, 'pade_extice', nband, nsizereg, ncoeff_ext, nrghice)
    pade_ssaice  = read_field(ncid, 'pade_ssaice', nband, nsizereg, ncoeff_ssa_g, nrghice)
    pade_asyice  = read_field(ncid, 'pade_asyice', nband, nsizereg, ncoeff_ssa_g, nrghice)

    ! Allocate cloud property Pade coefficient particle size boundary input arrays
    allocate(pade_sizreg_extliq(nbound), &
             pade_sizreg_ssaliq(nbound), &
             pade_sizreg_asyliq(nbound), &
             pade_sizreg_extice(nbound), &
             pade_sizreg_ssaice(nbound), &
             pade_sizreg_asyice(nbound))

    pade_sizreg_extliq = read_field(ncid, 'pade_sizreg_extliq', nbound)
    pade_sizreg_ssaliq = read_field(ncid, 'pade_sizreg_ssaliq', nbound)
    pade_sizreg_asyliq = read_field(ncid, 'pade_sizreg_asyliq', nbound)
    pade_sizreg_extice = read_field(ncid, 'pade_sizreg_extice', nbound)
    pade_sizreg_ssaice = read_field(ncid, 'pade_sizreg_ssaice', nbound)
    pade_sizreg_asyice = read_field(ncid, 'pade_sizreg_asyice', nbound)

    ncid = nf90_close(ncid)

    call stop_on_err(cloud_spec%load(band_lims_wvn, &
                                     pade_extliq, pade_ssaliq, pade_asyliq, &
                                     pade_extice, pade_ssaice, pade_asyice, &
                                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice))
  end subroutine load_cld_padecoeff

  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_cldpp(filename, cldfrac, clwp, ciwp, rei, rel, icergh)
    character(len=*),                      intent(in   ) :: filename
    real(wp), dimension(:,:), allocatable, intent(inout) :: cldfrac    ! cloud fraction
    real(wp), dimension(:,:), allocatable, intent(inout) :: ciwp       ! cloud ice water path
    real(wp), dimension(:,:), allocatable, intent(inout) :: clwp       ! cloud liquid water path
    real(wp), dimension(:,:), allocatable, intent(inout) :: rei        ! cloud ice particle effective size (microns)
    real(wp), dimension(:,:), allocatable, intent(inout) :: rel        ! cloud liquid particle effective radius (microns)
    integer,                               intent(  out) :: icergh
    ! -----------------
    integer :: ncid
    integer :: ncol, nlay
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("read_cldpp(): can't open file " // trim(filename))

    ncol = get_dim_size(ncid,'col')
    nlay = get_dim_size(ncid,'lay')
    ! Allocate cloud physical property arrays
    allocate(cldfrac(ncol,nlay), &
             clwp(ncol,nlay),    &
             ciwp(ncol,nlay),    &
             rel(ncol,nlay),     &
             rei(ncol,nlay))

    icergh     = read_field(ncid, 'icergh')
    cldfrac    = read_field(ncid, 'cloudfrac', ncol, nlay)
    clwp       = read_field(ncid, 'clwp',      ncol, nlay)
    ciwp       = read_field(ncid, 'ciwp',      ncol, nlay)
    rel        = read_field(ncid, 'rel',       ncol, nlay)
    rei        = read_field(ncid, 'rei',       ncol, nlay)

    ncid = nf90_close(ncid)
  end subroutine read_cldpp
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write cloud optical properties (nbnd)
  !
  subroutine write_cldop(filename, nbnd, cloud_optical_props)
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: nbnd
    class(ty_optical_props_arry), intent(in) :: cloud_optical_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nmom
    ! -------------------
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_cldop: can't open file " // trim(filename))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol   = get_dim_size(ncid, 'col')
    nlay   = get_dim_size(ncid, 'lay')
    nmom   = get_dim_size(ncid, 'mom')

    select type(cloud_optical_props)
      type is (ty_optical_props_1scl)
        call create_var(ncid,    "taucld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "taucld",  cloud_optical_props%tau))
      type is (ty_optical_props_2str)
        call create_var(ncid,    "taucld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "ssacld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "asycld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "taucld",  cloud_optical_props%tau))
        call stop_on_err(write_field(ncid, "ssacld",  cloud_optical_props%ssa))
        call stop_on_err(write_field(ncid, "asycld",  cloud_optical_props%g  ))
      type is (ty_optical_props_nstr)
        call create_var(ncid,    "taucld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "ssacld",  ["col ", "lay ", "band"], [ncol, nlay, nbnd])
        call create_var(ncid,    "pcld",    ["mom ", "col ", "lay ", "band"], [nmom, ncol, nlay, nbnd])
        call stop_on_err(write_field(ncid, "taucld",  cloud_optical_props%tau))
        call stop_on_err(write_field(ncid, "ssacld",  cloud_optical_props%ssa))
        call stop_on_err(write_field(ncid, "pcld",    cloud_optical_props%p  ))
    end select

    ncid = nf90_close(ncid)
  end subroutine write_cldop

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read cloud optical properties (nbnd) - for cloud sampling
  !
  !--------------------------------------------------------------------------------------------------------------------
!  subroutine read_cldop_nbnd(cloud_optical_props, opt_props, is_lw, filename, play, plev, tlay, &
  subroutine read_cldop_nbnd(cloud_optical_props, opt_props, filename, play, plev, tlay, &
                             ngpt, icldovrlp)

    class(ty_optical_props_arry),    intent(inout) :: cloud_optical_props
    class(ty_optical_props),               intent(inout) :: opt_props
!    logical,                               intent(in   ) :: is_lw
    character(len=*),                      intent(in   ) :: filename
    real(wp), dimension(:,:), allocatable, intent(  out) :: play, plev, tlay
    integer,                               intent(  out) :: ngpt
    integer,                               intent(  out) :: icldovrlp

    ! -----------------
    integer :: ncid
    integer :: ncol, nlay, nlev, nmom, nbnd, npai
    logical :: is_initialized

    integer, dimension(:,:), allocatable                 :: band_lims_gpt
    real(wp), dimension(:,:), allocatable                :: band_lims_wvn

    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("read_cldop_nbnd(): can't open file " // trim(filename))

    ncol = get_dim_size(ncid,'col')
    nlay = get_dim_size(ncid,'lay')
    nlev = get_dim_size(ncid,'lev')
    nmom = get_dim_size(ncid,'mom')
    npai = get_dim_size(ncid,'pair')

    icldovrlp = read_field(ncid, 'icldovrlp')

    allocate(play(ncol,nlay), plev(ncol,nlev), tlay(ncol,nlay))
    play = read_field(ncid, 'p_lay', ncol, nlay)
    plev = read_field(ncid, 'p_lev', ncol, nlev)
    tlay = read_field(ncid, 't_lay', ncol, nlay)

    if (is_lw(trim(filename))) then
       nbnd = get_dim_size(ncid, 'band')
       ngpt = get_dim_size(ncid, 'gpt')
       select type(cloud_optical_props)
         type is (ty_optical_props_1scl)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
         type is (ty_optical_props_2str)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
            cloud_optical_props%ssa    = read_field(ncid, 'ssacld', ncol, nlay, nbnd)
            cloud_optical_props%g      = read_field(ncid, 'asycld', ncol, nlay, nbnd)
         type is (ty_optical_props_nstr)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
            cloud_optical_props%ssa    = read_field(ncid, 'ssacld', ncol, nlay, nbnd)
            cloud_optical_props%p      = read_field(ncid, 'pcld', nmom, ncol, nlay, nbnd)
       end select
    end if

    if (.not. is_lw(trim(filename))) then
       nbnd = get_dim_size(ncid, 'band')
       ngpt = get_dim_size(ncid, 'gpt')
       select type(cloud_optical_props)
         type is (ty_optical_props_1scl)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
         type is (ty_optical_props_2str)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
            cloud_optical_props%ssa    = read_field(ncid, 'ssacld', ncol, nlay, nbnd)
            cloud_optical_props%g      = read_field(ncid, 'asycld', ncol, nlay, nbnd)
         type is (ty_optical_props_nstr)
            cloud_optical_props%tau    = read_field(ncid, 'taucld', ncol, nlay, nbnd)
            cloud_optical_props%ssa    = read_field(ncid, 'ssacld', ncol, nlay, nbnd)
            cloud_optical_props%p      = read_field(ncid, 'pcld', nmom, ncol, nlay, nbnd)
       end select
    end if

    ! Check for ty_optical_props initialization
    ! If false, initialize gpt2band here for use in cloud sampling
    if (.not. opt_props%is_initialized()) then
       allocate(band_lims_wvn(npai,nbnd))
       allocate(band_lims_gpt(npai,nbnd))
       band_lims_wvn = read_field(ncid, 'band_lims_wvn', npai, nbnd)
       band_lims_gpt = read_field(ncid, 'band_lims_gpt', npai, nbnd)
       call stop_on_err(opt_props%init(band_lims_wvn))
    end if

    ncid = nf90_close(ncid)

  end subroutine read_cldop_nbnd

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write cloud optical properties (ngpt)
  !
!  subroutine write_cldop_ngpt(filename, is_lw, cloud_optical_props)
  subroutine write_cldop_ngpt(filename, cloud_optical_props)
    character(len=*),                   intent(in) :: filename
!    logical,                            intent(in) :: is_lw
    class(ty_optical_props_arry), intent(in) :: cloud_optical_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngptlw, ngptsw, nmom
    ! -------------------
    if(nf90_open(trim(filename), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_cldop_ngpt: can't open file " // trim(filename))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol   = get_dim_size(ncid, 'col')
    nlay   = get_dim_size(ncid, 'lay')
    nmom   = get_dim_size(ncid, 'mom')

    if (is_lw(trim(filename))) then
       ngptlw = get_dim_size(ncid, 'gpt')
       select type(cloud_optical_props)
         type is (ty_optical_props_1scl)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
         type is (ty_optical_props_2str)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "ssacloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "asycloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssacloud",  cloud_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "asycloud",  cloud_optical_props%g ))
         type is (ty_optical_props_nstr)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "ssacloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptlw])
           call create_var(ncid,    "pcloud",    ["mom", "col", "lay", "gpt"], [nmom, ncol, nlay, ngptlw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssacloud",  cloud_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "pcloud",    cloud_optical_props%p ))
       end select
   endif

    if (.not. is_lw(trim(filename))) then
       ngptsw = get_dim_size(ncid, 'gpt')
       select type(cloud_optical_props)
         type is (ty_optical_props_1scl)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
         type is (ty_optical_props_2str)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "ssacloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "asycloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssacloud",  cloud_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "asycloud",  cloud_optical_props%g ))
         type is (ty_optical_props_nstr)
           call create_var(ncid,    "taucloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "ssacloud",  ["col", "lay", "gpt"], [ncol, nlay, ngptsw])
           call create_var(ncid,    "pcloud",    ["mom", "col", "lay", "gpt"], [nmom, ncol, nlay, ngptsw])
           call stop_on_err(write_field(ncid, "taucloud",  cloud_optical_props%tau ))
           call stop_on_err(write_field(ncid, "ssacloud",  cloud_optical_props%ssa ))
           call stop_on_err(write_field(ncid, "pcloud",    cloud_optical_props%p ))
       end select
    endif

    ncid = nf90_close(ncid)
  end subroutine write_cldop_ngpt

  !------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do SW calculations ?
  !
  function is_sw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_sw

    integer :: ncid

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("is_sw: can't find file " // trim(fileName))

    is_sw = var_exists(ncid, 'solar_zenith_angle')
    ncid = nf90_close(ncid)
  end function is_sw

  !------------------------------------------------------------------------------------------------------
  function is_lw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw

    is_lw = .not. is_sw(fileName)
  end function is_lw

  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        stop
      end if
    end subroutine

end module
