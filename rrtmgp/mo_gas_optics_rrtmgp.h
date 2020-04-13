
#pragma once

// This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
//
// Contacts: Robert Pincus and Eli Mlawer
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
// -------------------------------------------------------------------------------------------------
//
// Class for computing spectrally-resolved gas optical properties and source functions
//   given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
//   The class must be initialized with data (provided as a netCDF file) before being used.
//
// Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
//   external stellar radiation (shortwave radiation in the Earth's atmosphere).
//   The variant is chosen based on what information is supplied during initialization.
//   (It might make more sense to define two sub-classes)
//
// -------------------------------------------------------------------------------------------------

class GasOpticsRRTMGP : public OpticalProps {
public:
  // RRTMGP computes absorption in each band arising from
  //   two major species in each band, which are combined to make
  //     a relative mixing ratio eta and a total column amount (col_mix)
  //   contributions from zero or more minor species whose concentrations
  //     may be scaled by other components of the atmosphere
  // Absorption coefficients are interpolated from tables on a pressure/temperature/(eta) grid
  // Interpolation variables: Temperature and pressure grids
  real1d press_ref;
  real1d press_ref_log;
  real1d temp_ref;
  
  // Derived and stored for convenience:
  //   Min and max for temperature and pressure intepolation grids
  //   difference in ln pressure between consecutive reference levels
  //   log of reference pressure separating the lower and upper atmosphere
  real press_ref_min;
  real press_ref_max;
  real temp_ref_min;
  real temp_ref_max;
  real press_ref_log_delta;
  real temp_ref_delta;
  real press_ref_trop_log;
  
  // Major absorbers ("key species")
  //   Each unique set of major species is called a flavor.
  // Names  and reference volume mixing ratios of major gases
  string1d gas_names;  // gas names
  real3d vmr_ref;      // vmr_ref(lower or upper atmosphere, gas, temp)
  
  // Which two gases are in each flavor? By index
  int2d flavor;        // major species pair; (2,nflav)
  
  // Which flavor for each g-point? One each for lower, upper atmosphere
  int2d gpoint_flavor; // flavor = gpoint_flavor(2, g-point)
  
  // Major gas absorption coefficients
  real4d kmajor;       //  kmajor(g-point,eta,pressure,temperature)
  
  // Minor species, independently for upper and lower atmospheres
  //   Array extents in the n_minor dimension will differ between upper and lower atmospheres
  //   Each contribution has starting and ending g-points
  int2d minor_limits_gpt_lower;
  int2d minor_limits_gpt_upper;
  
  // Minor gas contributions might be scaled by other gas amounts; if so we need to know
  //   the total density and whether the contribution is scaled by the partner gas
  //   or its complement (i.e. all other gases)
  // Water vapor self- and foreign continua work like this, as do
  //   all collision-induced abosption pairs
  bool1d minor_scales_with_density_lower;
  bool1d minor_scales_with_density_upper;
  bool1d scale_by_complement_lower;
  bool1d scale_by_complement_upper;
  int1d idx_minor_lower;
  int1d idx_minor_upper;
  int1d idx_minor_scaling_lower;
  int1d idx_minor_scaling_upper;
  
  // Index into table of absorption coefficients
  int1d kminor_start_lower;
  int1d kminor_start_upper;
  
  // The absorption coefficients themselves
  real3d kminor_lower; // kminor_lower(n_minor,eta,temperature)
  real3d kminor_upper; // kminor_upper(n_minor,eta,temperature)

  // Rayleigh scattering coefficients
  real4d krayl; // krayl(g-point,eta,temperature,upper/lower atmosphere)

  // Planck function spectral mapping
  //   Allocated only when gas optics object is internal-source
  real4d planck_frac;   // stored fraction of Planck irradiance in band for given g-point
                        // planck_frac(g-point, eta, pressure, temperature)
  real2d totplnk;       // integrated Planck irradiance by band; (Planck temperatures,band)
  real   totplnk_delta; // temperature steps in totplnk
  
  // Solar source function spectral mapping
  //   Allocated only when gas optics object is external-source
  real1d solar_src; // incoming solar irradiance(g-point)

  // Ancillary
  // Index into %gas_names -- is this a key species in any band?
  bool1d is_key;


  // Everything except GasConcs is on the host by default, and the available_gases.gas_name is on the host as well
  // Things will be copied to the GPU outside of this routine and stored into class variables
  void reduce_minor_arrays(GasConcs   const &available_gases, 
                           string1d   const &gas_names, 
                           string1d   const &gas_minor,
                           string1d   const &identifier_minor,
                           realHost3d const &kminor_atm, 
                           string1d   const &minor_gases_atm, 
                           intHost2d  const &minor_limits_gpt_atm, 
                           boolHost1d const &minor_scales_with_density_atm, 
                           string1d   const &scaling_gas_atm, 
                           boolHost1d const &scale_by_complement_atm, 
                           intHost1d  const &kminor_start_atm, 
                           realHost3d       &kminor_atm_red, 
                           string1d         &minor_gases_atm_red, 
                           intHost2d        &minor_limits_gpt_atm_red, 
                           boolHost1d       &minor_scales_with_density_atm_red, 
                           string1d         &scaling_gas_atm_red, 
                           boolHost1d       &scale_by_complement_atm_red, 
                           intHost1d        &kminor_start_atm_red) {

    int nm = size(minor_gases_atm,1);  // Size of the larger list of minor gases
    int tot_g = 0;                     // TODO: Not documented in the original code
    int red_nm = 0;                    // Reduced number of minor gasses (only the ones we need)
    boolHost1d gas_is_present("gas_is_present",1,nm);   // Determines whether a gas in the list is needed
    // Determine the gasses needed
    for (int i=1; i <= nm; i++) {
      int idx_mnr = string_loc_in_array(minor_gases_atm(i), identifier_minor);
      gas_is_present(i) = string_in_array(gas_minor(idx_mnr),available_gases.gas_name);
      if (gas_is_present(i)) {
        tot_g = tot_g + (minor_limits_gpt_atm(2,i)-minor_limits_gpt_atm(1,i)+1);
        red_nm = red_nm + 1;
      }
    }

    // Allocate reduced arrays
    minor_gases_atm_red               = string1d  ("minor_gases_atm_red              "  ,red_nm);
    minor_scales_with_density_atm_red = boolHost1d("minor_scales_with_density_atm_red"  ,red_nm);
    scaling_gas_atm_red               = string1d  ("scaling_gas_atm_red              "  ,red_nm);
    scale_by_complement_atm_red       = boolHost1d("scale_by_complement_atm_red      "  ,red_nm);
    kminor_start_atm_red              = intHost1d ("kminor_start_atm_red             "  ,red_nm);
    minor_limits_gpt_atm_red          = intHost2d ("minor_limits_gpt_atm_red         ",2,red_nm);
    kminor_atm_red                    = realHost3d("kminor_atm_red                   ",tot_g , size(kminor_atm,2) , size(kminor_atm,3));

    if (red_nm == nm) {
      // If the gasses listed exactly matches the gasses needed, just copy it
      minor_gases_atm              .deep_copy_to(minor_gases_atm_red              );
      scaling_gas_atm              .deep_copy_to(scaling_gas_atm_red              );
      kminor_atm                   .deep_copy_to(kminor_atm_red                   );
      minor_limits_gpt_atm         .deep_copy_to(minor_limits_gpt_atm_red         );
      minor_scales_with_density_atm.deep_copy_to(minor_scales_with_density_atm_red);
      scale_by_complement_atm      .deep_copy_to(scale_by_complement_atm_red      );
      kminor_start_atm             .deep_copy_to(kminor_start_atm_red             );
    } else {
      // Otherwise, pack into reduced arrays
      int slot = 1;
      for (int i=1; i <= nm; i++) {
        if (gas_is_present(i)) {
          minor_gases_atm_red              (slot) = minor_gases_atm              (i);
          scaling_gas_atm_red              (slot) = scaling_gas_atm              (i);
          minor_scales_with_density_atm_red(slot) = minor_scales_with_density_atm(i);
          scale_by_complement_atm_red      (slot) = scale_by_complement_atm      (i);
          kminor_start_atm_red             (slot) = kminor_start_atm             (i);
          slot++;
        }
      }

      slot = 0;
      int n_elim = 0;
      for (int i=1 ; i <= nm ; i++) {
        int ng = minor_limits_gpt_atm(2,i)-minor_limits_gpt_atm(1,i)+1;
        if (gas_is_present(i)) {
          slot = slot + 1;
          minor_limits_gpt_atm_red   (1,slot) = minor_limits_gpt_atm(1,i);
          minor_limits_gpt_atm_red   (2,slot) = minor_limits_gpt_atm(2,i);
          kminor_start_atm_red         (slot) = kminor_start_atm(i)-n_elim;
          for (int l=1 ; l <= size(kminor_atm,3) ; l++ ) {
            for (int k=1 ; k <= size(kminor_atm,2) ; k++ ) {
              for (int j=1 ; j <= ng ; j++) {
                kminor_atm_red(kminor_start_atm_red(slot)+j-1,k,l) = kminor_atm(kminor_start_atm(i)+j-1,k,l);
              }
            }
          }
        } else {
          n_elim = n_elim + ng;
        }
      }
    }
  }



  // create index list for extracting col_gas needed for minor gas optical depth calculations
  void create_idx_minor(string1d const &gas_names, string1d const &gas_minor, string1d const &identifier_minor,
                        string1d const &minor_gases_atm, intHost1d &idx_minor_atm) {
    idx_minor_atm = intHost1d("idx_minor_atm",size(minor_gases_atm,1));
    for (int imnr=1 ; imnr <= size(minor_gases_atm,1) ; imnr++) {
      // Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
      int idx_mnr     = string_loc_in_array(minor_gases_atm(imnr), identifier_minor);
      // Find name of gas associated with minor species identifier (e.g. h2o)
      idx_minor_atm(imnr) = string_loc_in_array(gas_minor(idx_mnr), gas_names);
    }
  }



  // create index for special treatment in density scaling of minor gases
  void create_idx_minor_scaling(string1d const &gas_names, string1d const &scaling_gas_atm,
                                intHost1d &idx_minor_scaling_atm) {
    idx_minor_scaling_atm = intHost1d("idx_minor_scaling_atm",size(scaling_gas_atm,1));
    do imnr = 1, size(scaling_gas_atm,dim=1) ! loop over minor absorbers in each band
    for (int imnr=1 ; imnr <= size(scaling_gas_atm,1) ; imnr++) {
      // This will be -1 if there's no interacting gas
      idx_minor_scaling_atm(imnr) = string_loc_in_array(scaling_gas_atm(imnr), gas_names);
    }
  }



  void create_key_species_reduce(string1d const &gas_names, string1d const &gas_names_red, intHost3d const &key_species,
                                 intHost3d &key_species_red, boolHost1d &key_species_present_init) {
    int np = size(key_species,1);
    int na = size(key_species,2);
    int nt = size(key_species,3);
    key_species_red = intHost3d("key_species_red",np,na,nt);
    key_species_present_init = boolHost1d("key_species_present_init",size(gas_names,1));

    for (int i=1 ; i <= size(gase_names,1) ; i++) {
      key_species_present_init(i) = true;
    }

    for (int ip=1 ; ip <= np ; ip++) {
      for (int ia=1 ; ia <= na ; ia++) {
        for (int it=1 ; it <= nt ; it++) {
          if (key_species(ip,ia,it) != 0) {
            key_species_red(ip,ia,it) = string_loc_in_array(gas_names(key_species(ip,ia,it)),gas_names_red);
            if (key_species_red(ip,ia,it) == -1) {
              key_species_present_init(key_species(ip,ia,it)) = false;
            }
          } else {
            key_species_red(ip,ia,it) = key_species(ip,ia,it);
          }
        }
      }
    }
  }



  // Create flavor list
  // An unordered array of extent (2,:) containing all possible pairs of key species used in either upper or lower atmos
  void create_flavor(intHost3d const &key_species, intHost2d &flavor) {
    // prepare list of key_species
    int i = 1;
    intHost2d key_species_list("key_species_list",2,size(key_species,3)*2);
    for (int ibnd=1 ; ibnd <= size(key_species,3) ; ibnd++) {
      for (int iatm=1 ; iatm <= size(key_species,1) ; iatm++) {
        key_species_list(1,i) = key_species(1,iatm,ibnd);
        key_species_list(2,i) = key_species(2,iatm,ibnd);
        i = i + 1;
      }
    }
    // rewrite single key_species pairs
    for (int i=1 ; i = size(key_species_list,2) ; i++) {
      if (key_species_list(1,i) == 0 && key_species_list(2,i) == 0) {
        key_species_list(1,i) = 2;
        key_species_list(2,i) = 2;
      }
    }
    // count unique key species pairs
    int iflavor = 0;
    for (int i=1; i <= size(key_species_list,2) ; i++) {
      // Loop through previous pairs. Only increment iflavor if we haven't seen this pair before
      bool unique = true;
      for (int j=1; j <= i-1 ; j++) {
        if ( key_species_list(1,j) == key_species_list(1,i) && key_species_list(2,j) == key_species_list(2,i) ) {
          unique = false;
        }
      }
      if (unique) { iflavor = iflavor + 1; }
    }
    // fill flavors
    flavor = intHost2d("flavor",2,iflavor);
    iflavor = 0;
    for (int i=1 ; i <= size(key_species_list,2) ; i++) {
      bool unique = true;
      for (int j=1; j <= i-1 ; j++) {
        if ( key_species_list(1,j) == key_species_list(1,i) && key_species_list(2,j) == key_species_list(2,i) ) {
          unique = false;
        }
      }
      if (unique) {
        iflavor = iflavor + 1;
        flavor(1,iflavor) = key_species_list(1,i);
        flavor(2,iflavor) = key_species_list(2,i);
      }
    }
  }



  // returns flavor index; -1 if not found
  pure function key_species_pair2flavor(flavor, key_species_pair)
    integer :: key_species_pair2flavor
    integer, dimension(:,:), intent(in) :: flavor
    integer, dimension(2), intent(in) :: key_species_pair
    integer :: iflav
    do iflav=1,size(flavor,2)
      if (all(key_species_pair(:).eq.flavor(:,iflav))) then
        key_species_pair2flavor = iflav
        return
      end if
    end do
    key_species_pair2flavor = -1
  end function key_species_pair2flavor



  // create gpoint_flavor list: a map pointing from each g-point to the corresponding entry in the "flavor list"
  void create_gpoint_flavor(intHost3d const &key_species, intHost1d const &gpt2band, intHost2d const &flavor,
                            intHost2d &gpoint_flavor) {
    int ngpt = size(gpt2band,1);
    gpoint_flavor = intHost2d("gpoint_flavor",2,ngpt);
    for (int igpt=1 ; igpt <= ngpt ; igpt++) {
      for (int iatm = 1 ; iatm <= 2 ; iatm++) {
        int key_species_pair2flavor = -1;
        for (int iflav=1 ; iflav <= size(flavor,2) ; iflav++) {
          int key_species1 = key_species(1,iatm,gpt2band(igpt));
          int key_species2 = key_species(2,iatm,gpt2band(igpt));
          if (key_species1 == 0 && key_species2 == 0) {
            key_species1 = 2;
            key_species2 = 2;
          }
          if ( flavor(1,iflav) == key_species1 && flavor(2,iflav) == key_species2 ) {
            key_species_pair2flavor = iflav;
          }
        }
        gpoint_flavor(iatm,igpt) = key_speciespair2flavor;
      }
    }
  }



  // Initialize absorption coefficient arrays,
  //   including Rayleigh scattering tables if provided (allocated)
  void init_abs_coeffs(GasConcs   const &available_gases, 
                       string1d   const &gas_names,
                       intHost3d  const &key_species,    
                       intHost2d  const &band2gpt,
                       realHost2d const &band_lims_wavenum, 
                       realHost1d const &press_ref,
                       realHost1d const &temp_ref,       
                       real              press_ref_trop,
                       real              temp_ref_p,
                       real              temp_ref_t, 
                       realHost3d const &vmr_ref,                   
                       realHost4d const &kmajor,
                       realHost3d const &kminor_lower,
                       realHost3d const &kminor_upper, 
                       string1d   const &gas_minor,
                       string1d   const &identifier_minor,
                       string1d   const &minor_gases_lower,
                       string1d   const &minor_gases_upper, 
                       intHost2d  const &minor_limits_gpt_lower, 
                       intHost2d  const &minor_limits_gpt_upper, 
                       boolHost1d const &minor_scales_with_density_lower, 
                       boolHost1d const &minor_scales_with_density_upper, 
                       string1d   const &scaling_gas_lower,
                       string1d   const &scaling_gas_upper, 
                       boolHost1d const &scale_by_complement_lower, 
                       boolHost1d const &scale_by_complement_upper, 
                       intHost1d  const &kminor_start_lower, 
                       intHost1d  const &kminor_start_upper, 
                       realHost3d const &rayl_lower,
                       realHost3d const &rayl_upper) {

    OpticalProps::init(band_lims_wavenum, band2gpt);
    // Which gases known to the gas optics are present in the host model (available_gases)?
    int ngas = size(gas_names,1);
    boolHost1d gas_is_present("gas_is_present",ngas);
    for (int i=1 ; i <= ngas ; i++) {
      gas_is_present(i) = string_in_array(gas_names(i), available_gases.gas_name);
    }
    // Now the number of gases is the union of those known to the k-distribution and provided by the host model
    ngas = count(gas_is_present);
    
    // Initialize the gas optics object, keeping only those gases known to the gas optics and also present in the host model
    this->gas_names = pack(gas_names,gas_is_present);

    realHost3d vmr_ref_red("vmr_ref_red",size(vmr_ref,1),{0,ngas}, size(vmr_ref,3));
    // Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    for (int k=1 ; k <= size(vmr_ref,3) ; k++) {
      for (int j=1 ; j <= size(vmr_ref,1) ; j++) {
        vmr_ref_red(j,0,k) = vmr_ref(j,1,k);
      }
    }
    do i = 1, ngas
    for (int i=1 ; i <= ngas ; i++) {
      int idx = string_loc_in_array(this%gas_names(i), gas_names);
      for (int k=1 ; k <= size(vmr_ref,3) ; k++) {
        for (int j=1 ; j <= size(vmr_ref,1) ; j++) {
          vmr_ref_red(j,i,k) = vmr_ref(j,idx+1,k);
        }
      }
    }
    // Allocate class copy, and deep copy to the class data member
    this->vmr_ref = real3d("vmr_ref",size(vmr_ref,1),{0,ngas}, size(vmr_ref,3));
    vmr_ref_red.deep_copy_to(this->vmr_ref);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // REDUCE MINOR ARRAYS SO VARIABLES ONLY CONTAIN MINOR GASES THAT ARE AVAILABLE
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LOWER MINOR GASSES
    string1d minor_gases_lower_red;
    string1d scaling_gas_lower_red;
    realHost3d kminor_lower_red;
    intHost2d  minor_limits_gpt_lower_red;
    boolHost1d minor_scales_with_density_lower_red;
    boolHost1d scale_by_complement_lower_red;
    intHost1d  kminor_start_lower_red;

    reduce_minor_arrays(available_gases, gas_names, gas_minor, identifier_minor, kminor_lower, minor_gases_lower, 
                        minor_limits_gpt_lower, minor_scales_with_density_lower, scaling_gas_lower, 
                        scale_by_complement_lower, kminor_start_lower, kminor_lower_red, minor_gases_lower_red, 
                        minor_limits_gpt_lower_red, minor_scales_with_density_lower_red, scaling_gas_lower_red, 
                        scale_by_complement_lower_red, kminor_start_lower_red);

    // Copy Host temps to class device data members
    kminor_lower_red                   .deep_copy_to(this->kminor_lower                   );
    minor_limits_gpt_lower_red         .deep_copy_to(this->minor_limits_gpt_lower         );
    minor_scales_with_density_lower_red.deep_copy_to(this->minor_scales_with_density_lower);
    scale_by_complement_lower_red      .deep_copy_to(this->scale_by_complement_lower      );
    kminor_start_lower_red             .deep_copy_to(this->kminor_start_lower             );

    // UPPER MINOR GASSES
    string1d minor_gases_upper_red;
    string1d scaling_gas_upper_red;
    realHost3d kminor_upper_red;
    intHost2d  minor_limits_gpt_upper_red;
    boolHost1d minor_scales_with_density_upper_red;
    boolHost1d scale_by_complement_upper_red;
    intHost1d  kminor_start_upper_red;

    call reduce_minor_arrays(available_gases, 
                             gas_names, 
                             gas_minor,identifier_minor,
                             kminor_upper, 
                             minor_gases_upper, 
                             minor_limits_gpt_upper, 
                             minor_scales_with_density_upper, 
                             scaling_gas_upper, 
                             scale_by_complement_upper, 
                             kminor_start_upper, 
                             kminor_upper_red, 
                             minor_gases_upper_red, 
                             minor_limits_gpt_upper_red, 
                             minor_scales_with_density_upper_red, 
                             scaling_gas_upper_red, 
                             scale_by_complement_upper_red, 
                             kminor_start_upper_red)

    // Copy Host temps to class device data members
    kminor_upper_red                   .deep_copy_to(this->kminor_upper                   );
    minor_limits_gpt_upper_red         .deep_copy_to(this->minor_limits_gpt_upper         );
    minor_scales_with_density_upper_red.deep_copy_to(this->minor_scales_with_density_upper);
    scale_by_complement_upper_red      .deep_copy_to(this->scale_by_complement_upper      );
    kminor_start_upper_red             .deep_copy_to(this->kminor_start_upper             );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // HANDLE ARRAYS NOT REDUCED BY THE PRESENCE, OR LACK THEREOF, OF A GAS
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    this->press_ref = real1d("press_ref",size(press_ref,1));
    this->temp_ref  = real1d("temp_ref" ,size(temp_ref ,1));
    this->kmajor    = real4d("kmajor"   ,size(kmajor,1),size(kmajor,2),size(kmajor,3),size(kmajor,4));
    press_ref.deep_copy_to(this->press_ref);
    temp_ref .deep_copy_to(this->temp_ref );
    kmajor   .deep_copy_to(this->kmajor   );

    // Process rayl_lower and rayl_upper into a combined this->krayl
    if (allocated(rayl_lower) != allocated(rayl_upper)) {
      stoprun("rayl_lower and rayl_upper must have the same allocation status");
    }
    if (allocated(rayl_lower)) {
      realHost4d krayltmp("krayltmp",size(rayl_lower,1),size(rayl_lower,2),size(rayl_lower,3),2);
      for (int k=1 ; k <= size(rayl_lower,3) ; k++ ) {
        for (int j=1 ; j <= size(rayl_lower,2) ; j++ ) {
          for (int i=1 ; i <= size(rayl_lower,1) ; i++ ) {
            krayltmp(i,j,k,1) = rayl_lower(i,j,k);
            krayltmp(i,j,k,2) = rayl_upper(i,j,k);
          }
        }
      }
      krayltmp.deep_copy_to(this->krayl);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // POST PROCESSING
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // creates log reference pressure
    this->press_ref_log = real1d("press_ref_log",size(this->press_ref,1));
    auto &press_ref_loc     = this->press_ref;
    auto &press_ref_log_loc = this->press_ref_log;
    // Running a kernel because it's more convenient in this case
    parallel_for_cpu_serial( Bounds<1>( size(this->press_ref,1) ) , YAKL_LAMBDA (int i) {
      press_ref_log_loc(i) = log(press_ref_loc(i));
    });

    // log scale of reference pressure (this is a scalar, not an array)
    this->press_ref_trop_log = log(press_ref_trop);

    // Get index of gas (if present) for determining col_gas
    intHost1d idx_minor_lower_tmp;
    intHost1d idx_minor_upper_tmp;
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_lower_red, idx_minor_lower_tmp);
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_upper_red, idx_minor_upper_tmp);
    idx_minor_lower_tmp.deep_copy_to(this->idx_minor_lower);
    idx_minor_upper_tmp.deep_copy_to(this->idx_minor_upper);
    // Get index of gas (if present) that has special treatment in density scaling
    intHost1d idx_minor_scaling_lower_tmp;
    intHost1d idx_minor_scaling_upper_tmp;
    create_idx_minor_scaling(this%gas_names, scaling_gas_lower_red, idx_minor_scaling_lower_tmp);
    create_idx_minor_scaling(this%gas_names, scaling_gas_upper_red, idx_minor_scaling_upper_tmp);
    idx_minor_scaling_lower_tmp.deep_copy_to(this->idx_minor_scaling_lower);
    idx_minor_scaling_upper_tmp.deep_copy_to(this->idx_minor_scaling_upper);

    // create flavor list
    // Reduce (remap) key_species list; checks that all key gases are present in incoming
    boolHost1d key_species_present_init;
    intHost2d  key_species_red;
    create_key_species_reduce(gas_names, this->gas_names, key_species, key_species_red, key_species_present_init);
    // create flavor and gpoint_flavor lists
    intHost2d flavor_tmp;
    intHost2d gpoint_flavor_tmp;
    create_flavor       (key_species_red, flavor_tmp);
    create_gpoint_flavor(key_species_red, this->get_gpoint_bands(), flavor_tmp, gpoint_flavor_tmp);
    this->flavor        = int2d("flavor",size(       flavor_tmp,1),size(       flavor_tmp,2));
    this->gpoint_flavor = int2d("flavor",size(gpoint_flavor_tmp,1),size(gpoint_flavor_tmp,2));
    flavor_tmp       .deep_copy_to(this->flavor       );
    gpoint_flavor_tmp.deep_copy_to(this->gpoint_flavor);

    // minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    //   for T, high-to-low ordering for p
    this->temp_ref_min  = temp_ref (1);
    this->temp_ref_max  = temp_ref (size(temp_ref ,1));
    this->press_ref_min = press_ref(size(press_ref,1));
    this->press_ref_max = press_ref(1);

    // creates press_ref_log, temp_ref_delta
    this->press_ref_log_delta = (log(this->press_ref_min)-log(this->press_ref_max))/(size(press_ref,1)-1);
    this->temp_ref_delta      = (this->temp_ref_max-this->temp_ref_min)/(size(temp_ref,1)-1);

    // Which species are key in one or more bands?
    //   this->flavor is an index into this->gas_names
    this->is_key = bool1d("is_key",this->get_ngas());
    auto &is_key_loc = this->is_key;
    auto &flavor_loc = this->flavor;
    parallel_for_cpu_serial( Bounds<1>( this->get_ngas() ) , YAKL_LAMBDA (int i) {
      is_key_loc(i) = false;
    };
    // do j = 1, size(this%flavor, 2)
    //   do i = 1, size(this%flavor, 1) ! extents should be 2
    parallel_for_cpu_serial( Bounds<2>( size(this->flavor,2) , size(this->flavor,1) ) , YAKL_LAMBDA (int j, int i) {
      if (flavor_loc(i,j) /= 0) { is_key_loc(flavor_loc(i,j)) = true; }
    });
  }



  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the internal-sources object -- includes Plank functions and fractions
  void load(GasConcs   const &available_gases,
            string1d   const &gas_names,
            intHost3d  const &key_species,  
            intHost2d  const &band2gpt,
            realHost2d const &band_lims_wavenum,                    
            realHost1d const &press_ref,
            real              press_ref_trop,
            realHost1d const &temp_ref,            
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3d const &vmr_ref,                
            realHost4d const &kmajor,
            realHost3d const &kminor_lower,
            realHost3d const &kminor_upper,             
            string1d   const &gas_minor,
            string1d   const &identifier_minor,                     
            string1d   const &minor_gases_lower,
            string1d   const &minor_gases_upper,           
            intHost2d  const &minor_limits_gpt_lower,
            intHost2d  const &minor_limits_gpt_upper, 
            boolHost1d const &minor_scales_with_density_lower,                
            boolHost1d const &minor_scales_with_density_upper,                
            string1d   const &scaling_gas_lower,
            string1d   const &scaling_gas_upper,           
            boolHost1d const &scale_by_complement_lower,                      
            boolHost1d const &scale_by_complement_upper,                      
            intHost1d  const &kminor_start_lower,                             
            intHost1d  const &kminor_start_upper,                             
            realHost2d const &totplnk,
            realHost4d const &planck_frac,
            realHost3d const &rayl_lower,
            realHost3d const &rayl_upper) {

    init_abs_coeffs(available_gases, gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,       
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper, 
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower, 
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper, 
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper, 
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);

    // Planck function tables
    this->totplnk = real2d("totplnk",size(totplnk,1),size(totplnk,2));
    this->planck_frac = real4d("planck_frac",size(planck_frac,1),size(planck_frac,2),size(planck_frac,3),size(planck_frac,4));
    totplnk    .deep_copy_to(this->totplnk);
    planck_frac.deep_copy_to(this->planck_frac);
    // Temperature steps for Planck function interpolation
    //   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    //   Planck grid and the Planck grid is equally spaced
    this->totplnk_delta = (this->temp_ref_max - this->temp_ref_min) / (size(this->totplnk,1)-1);
  }



  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the external-sources object -- includes TOA source function table
  void load(GasConcs   const &available_gases,
            string1d   const &gas_names,
            intHost3d  const &key_species,
            intHost2d  const &band2gpt,
            realHost2d const &band_lims_wavenum,
            realHost1d const &press_ref,
            real              press_ref_trop,
            realHost1d const &temp_ref,
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3d const &vmr_ref,
            realHost4d const &kmajor,
            realHost3d const &kminor_lower,
            realHost3d const &kminor_upper,
            string1d   const &gas_minor,
            string1d   const &identifier_minor,
            string1d   const &minor_gases_lower,
            string1d   const &minor_gases_upper,
            intHost2d  const &minor_limits_gpt_lower,
            intHost2d  const &minor_limits_gpt_upper,
            boolHost1d const &minor_scales_with_density_lower,
            boolHost1d const &minor_scales_with_density_upper,
            string1d   const &scaling_gas_lower,
            string1d   const &scaling_gas_upper,
            boolHost1d const &scale_by_complement_lower,
            boolHost1d const &scale_by_complement_upper,
            intHost1d  const &kminor_start_lower,
            intHost1d  const &kminor_start_upper,
            realHost1d const &solar_src,
            realHost3d const &rayl_lower,
            realHost3d const &rayl_upper) {

    init_abs_coeffs(available_gases,  gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);
    
    // Solar source table init
    this->solar_src = real1d("solar_src",size(solar_src,1));
    solar_src.deep_copy_to(this->solar_src);
  }



  // Two functions to define array sizes needed by gas_optics()
  int get_ngas() { return size(this->gas_names,1); }



  // return the number of distinct major gas pairs in the spectral bands (referred to as
  // "flavors" - all bands have a flavor even if there is one or no major gas)
  int get_nflav() { return size(this->flavor,2); }



  string1d get_gases() { return this->gas_names; }



  // return the minimum pressure on the interpolation grids
  real get_press_min() { return this->press_ref_min; }



  // return the maximum pressure on the interpolation grids
  real get_press_max() { return this->press_ref_max; }



  // return the minimum temparature on the interpolation grids
  real get_temp_min() { return this->temp_ref_min; }



  // return the maximum temparature on the interpolation grids
  real get_temp_max() { return this->temp_ref_max; }



  int get_neta() { return size(this->kmajor,2); }



  // return the number of pressures in reference profile
  //   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
  int get_npres() { return size(this->kmajor,3)-1; }



  int get_ntemp() { return size(this->kmajor,4); }



  // return the number of temperatures for Planck function
  int get_nPlanckTemp() { return size(this->totplnk,1); }



  // Function to define names of key and minor gases to be used by gas_optics().
  // The final list gases includes those that are defined in gas_optics_specification
  // and are provided in ty_gas_concs.
  string1d get_minor_list(this, GasConcs const &gas_desc, int ngas, string1d const &names_spec) {
    // List of minor gases to be used in gas_optics()
    boolHost1d gas_is_present("gas_is_present",size(name_spec,1));
    for (int igas=1 ; igas <= this->get_ngas() ; igas++) {
      gas_is_present(igas) = string_in_array(names_spec(igas), gas_desc.gas_name)
    };
    return pack(this->gas_names, gas_is_present);
  }



  // return true if initialized for internal sources, false otherwise
  bool source_is_internal() { return allocated(this->totplnk) && allocated(this->planck_frac); }



  // return true if initialized for external sources, false otherwise
  bool source_is_external() { return allocated(this->solar_src); }



  // Ensure that every key gas required by the k-distribution is present in the gas concentration object
  void check_key_species_present(GasConcs const &gas_desc) {
    string1d key_gas_names = pack(this->gas_names, this->is_key);
    for (int igas=1 ; igas <= size(key_gas_names,1) ; igas++) {
      if (! string_in_array(key_gas_names(igas), gas_desc.gas_name)) {
        stoprun("gas required by k-distribution is not present in the GasConcs object");
      }
    }
  }



  // Compute gas optical depth and Planck source functions, given temperature, pressure, and composition
  void gas_optics_int(real2d const &play, real2d const &plev, real2d const &tlay, real1d const &tsfc,
                      GasConcs const &gas_desc, OpticalPropsArry &optical_props, SourceFuncLW &sources,
                      real2d const &col_dry=real2d(), real2d const &tlev=real2d()) {
    int ncol  = size(play,1);
    int nlay  = size(play,2);
    int ngpt  = this->get_ngpt();
    int nband = this->get_nband();
    // Interpolation coefficients for use in source function
    int2d  jtemp ("jtemp"                         ,size(play,1),size(play,2));
    int2d  jpress("jpress"                        ,size(play,1),size(play,2));
    bool2d tropo ("tropo"                         ,size(play,1),size(play,2));
    real6d fmajor("fmajor",2,2,2,this->get_nflav(),size(play,1),size(play,2));
    int4d  jeta  ("jeta"  ,2    ,this->get_nflav(),size(play,1),size(play,2));
    // Gas optics
    compute_gas_taus(ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress,
                     jeta, tropo, fmajor, col_dry);

    // External source -- check arrays sizes and values
    // input data sizes and values
    if (size(tsfc,1) != ncol) { stoprun("gas_optics(): array tsfc has wrong size"); }
    if (anyLT(tsfc,this->temp_ref_min) || anyGT(tsfc,this->temp_ref_max)) {
      stoprun("gas_optics(): array tsfc has values outside range");
    }

    if (allocated(tlev)) {
      if (size(tlev,1) != ncol || size(tlev,2) != nlay+1) { stoprun("gas_optics(): array tlev has wrong size"); }
      if (anyLT(tlev,this->temp_ref_min) || anyGT(tlev,this->temp_ref_max)) {
        stoprun("gas_optics(): array tlev has values outside range");
      }
    }

    // output extents
    if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
      stoprun("gas_optics%gas_optics: source function arrays inconsistently sized");
    }

    // Interpolate source function
    source(ncol, nlay, nband, ngpt, play, plev, tlay, tsfc, jtemp, jpress, jeta, tropo, fmajor, sources, tlev);
  }



  // Compute gas optical depth given temperature, pressure, and composition
  function gas_optics_ext(real2d const &play, real2d const &plev, real2d const &tlay, GasConcs const &gas_desc,   
                          OpticalPropsArry &optical_props, real2d &toa_src, real2d const &col_dry=real2d()) {
    int ncol  = size(play,dim=1)
    int nlay  = size(play,dim=2)
    int ngpt  = this%get_ngpt()
    int nband = this%get_nband()
    int ngas  = this%get_ngas()
    int nflav = get_nflav(this)
    
    // Interpolation coefficients for use in source function
    int2d  jtemp ("jtemp"                         ,size(play,1),size(play,2));
    int2d  jpress("jpress"                        ,size(play,1),size(play,2));
    bool2d tropo ("tropo"                         ,size(play,1),size(play,2));
    real6d fmajor("fmajor",2,2,2,this->get_nflav(),size(play,1),size(play,2));
    int4d  jeta  ("jeta  ",2    ,this->get_nflav(),size(play,1),size(play,2));
    // Gas optics
    compute_gas_taus(ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress, jeta,
                     tropo, fmajor, col_dry);

    // External source function is constant
    if (size(toa_src,1) != ncol || size(toa_src,2) != ngpt) { stoprun("gas_optics(): array toa_src has wrong size"); }

    auto &solar_source_loc = this->solar_src;
    parallel_for_cpu_serial( Bounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      toa_src(icol,igpt) = solar_src_loc(igpt);
    });
  }


};

