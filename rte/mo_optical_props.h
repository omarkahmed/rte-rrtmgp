
#pragma once

#include "const.h"

// Base class for optical properties
//   Describes the spectral discretization including the wavenumber limits
//   of each band (spectral region) and the mapping between g-points and bands
class OpticalProps {
public:
  int2d  band2gpt;       // (begin g-point, end g-point) = band2gpt(2,band)
  int1d  gpt2band;       // band = gpt2band(g-point)
  real2d band_lims_wvn;  // (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  std::string name;


  // Base class: Initialization
  //   Values are assumed to be defined in bands a mapping between bands and g-points is provided
  void init_base(real2d const &band_lims_wvn, real2d const &band_lims_gpt=real2d(), std::string name="") {
    int2d band_lims_gpt_lcl("band_lims_gpt_lcl",2,size(band_lims_wvn,2));

    // Error checking -- are the arrays the size we expect, contain positive values?
    if ( size(band_lims_wvn,1) != 2 ) { stoprun("OpticalProps::init(): band_lims_wvn 1st dim should be 2"); }
    if ( any_vals_less_than(band_lims_wvn,0.) ) { stoprun("OpticalProps::init(): band_lims_wvn has values <  0."); }
    if ( band_lims_gpt.initialized() ) {
      if ( size(band_lims_gpt,1) != 2 || size(band_limes_gpt,2) != size(band_lims_wvn,2) )  {
        stoprun("OpticalProps::init(): band_lims_gpt size inconsistent with band_lims_wvn");
      }
      if ( any_vals_less_than(band_lims_gpt,1.) ) {
        stoprun("OpticalProps::init(): band_lims_gpt has values < 1");
      }
      parallel_for_cpu_serial( Bounds<2>(size(band_lims_wvn,2),2) , YAKL_LAMBDA (int j, int i) {
        band_lims_gpt_lcl(i,j) = band_lims_gpt(i,j);
      });
    } else {
      // Assume that values are defined by band, one g-point per band
      parallel_for_cpu_serial( size(band_lims_wvn,2) , YAKL_LAMBDA (int iband) {
        band_lims_gpt_lcl(1,iband) = iband;
        band_lims_gpt_lcl(2,iband) = iband;
      });
    }
    // Assignment
    this->band2gpt      = band_lims_gpt_lcl
    this->band_lims_wvn = band_lims_wvn
    this->name          = name;

    // Make a map between g-points and bands
    //   Efficient only when g-point indexes start at 1 and are contiguous.
    this->gpt2band = int1d("gpt2band",maxval(band_lims_gpt_lcl));

    if(allocated(this%gpt2band)) deallocate(this%gpt2band)

    allocate(this%gpt2band(maxval(band_lims_gpt_lcl)))

    do iband=1,size(band_lims_gpt_lcl,dim=2)
      this%gpt2band(band_lims_gpt_lcl(1,iband):band_lims_gpt_lcl(2,iband)) = iband
    end do
  }



  // generic,   public  :: init => init_base, init_base_from_copy
  // procedure, private :: init_base
  // procedure, private :: init_base_from_copy
  // procedure, public  :: is_initialized => is_initialized_base
  // procedure, private :: is_initialized_base
  // procedure, public  :: finalize => finalize_base
  // procedure, private :: finalize_base
  // procedure, public  :: get_nband
  // procedure, public  :: get_ngpt
  // procedure, public  :: get_gpoint_bands
  // procedure, public  :: convert_band2gpt
  // procedure, public  :: convert_gpt2band
  // procedure, public  :: get_band_lims_gpoint
  // procedure, public  :: get_band_lims_wavenumber
  // procedure, public  :: get_band_lims_wavelength
  // procedure, public  :: bands_are_equal
  // procedure, public  :: gpoints_are_equal
  // procedure, public  :: expand
  // procedure, public  :: set_name
  // procedure, public  :: get_name
};




