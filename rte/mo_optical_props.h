
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
  void init( real2d &band_lims_wvn , int2d &band_lims_gpt=ind2d() , std::string name="" ) {
    int2d band_lims_gpt_lcl("band_lims_gpt_lcl",2,size(band_lims_wvn,2));
    if (size(band_lims_wvn,1) != 2) { stoprun("optical_props::init(): band_lims_wvn 1st dim should be 2"); }
    if (anyLT(band_lims_wvn,0._wp)) { stoprun("optical_props::init(): band_lims_wvn has values <  0."); }
    if (allocated(band_lims_gpt)) {
      if (size(band_lims_gpt,2) != size(band_lims_wvn,2)) {
        stoprun("optical_props::init(): band_lims_gpt size inconsistent with band_lims_wvn");
      }
      if (anyLT(band_lims_gpt,1._wp) ) { stoprun("optical_props%init(): band_lims_gpt has values < 1"); }
      // for (int j=1; j <= size(band_lims_gpt,2); j++) {
      //   for (int i=1; i <= size(band_lims_gpt,1); i++) {
      parallel_for_cpu_serial( Bounds<2>(size(band_lims_gpt,2),size(band_lims_gpt,1)) , YAKL_LAMBDA (int j, int i) {
        band_lims_gpt_lcl(i,j) = band_lims_gpt(i,j);
      });
    } else {
      // Assume that values are defined by band, one g-point per band
      // for (int iband = 1; iband <= size(band_lims_wvn, 2); iband++) {
      parallel_for_cpu_serial( Bounds<1>(size(band_lims_wvn, 2)) , YAKL_LAMBDA (int iband) {
        band_lims_gpt_lcl(2,iband) = iband;
        band_lims_gpt_lcl(1,iband) = iband;
      });
    }
    // Assignment
    this->band2gpt      = band_lims_gpt_lcl;
    this->band_lims_wvn = band_lims_wvn;
    this->name          = name;

    // Make a map between g-points and bands
    //   Efficient only when g-point indexes start at 1 and are contiguous.
    this->gpt2band = int1d("gpt2band",maxval(band_lims_gpt_lcl));
    // TODO: I didn't want to bother with race conditions at the moment, so it's an entirely serialized kernel for now
    parallel_for_cpu_serial( Bounds<1>(1) , YAKL_LAMBDA (int dummy) {
      for (int iband=1; iband <= size(band_lims_gpt_lcl,2); iband++) {
        for (int i=band_lims_gpt_lcl(1,iband); i <= band_lims_gpt_lcl(2,iband); i++) {
          this->gpt2band(i) = iband;
        }
      }
    });
  }


  void init(OptcalProps const &in) {
    if ( ! in.is_initialized() ) {
      stoprun("optical_props::init(): can't initialize based on un-initialized input");
    } else {
      this->init( in.get_band_lims_wavenumber() , in.get_band_lims_gpoint() );
    }
  }


  YAKL_INLINE bool is_initialized() { return allocated(this->band2gpt); }


  // Base class: finalize (deallocate memory)
  void finalize() {
    this->band2gpt      = int2d();
    this->gpt2band      = int1d();
    this->band_lims_wvn = real2d();
    this->name          = "";
  }


  // Number of bands
  YAKL_INLINE int get_nband() {
    if (this->is_initialized()) { return size(this->band2gpt,2); }
    return 0;
  }


  // Number of g-points
  YAKL_INLINE int get_ngpt()
    if (this->is_initialized()) { return maxval(this->band2gpt); }
    return 0;
  }


  // Bands for all the g-points at once;  dimension (ngpt)
  int1d get_gpoint_bands() { return gpt2band; }


  // First and last g-point of a specific band
  int1d convert_band2gpt(band) {
    int1d ret("band2gpt",2);
    if (this->is_initialized()) {
      ret(1) = this->band2gpt(1,band);
      ret(2) = this->band2gpt(2,band);
    } else {
      ret(1) = 0;
      ret(2) = 0;
    }
    return ret;
  }


  // Band associated with a specific g-point
  int convert_gpt2band(int gpt) {
    if (this->is_initialized()) { return this->gpt2band(gpt); }
    return 0;
  }


  // The first and last g-point of all bands at once;  dimension (2, nbands)
  int2d get_band_lims_gpoint() { return this%band2gpt; }


  // Lower and upper wavenumber of all bands
  // (upper and lower wavenumber by band) = band_lims_wvn(2,band)
  real2d get_band_lims_wavenumber() { return this->band_lims_wvn; }


  // Lower and upper wavelength of all bands
  real2d function get_band_lims_wavelength() {
    real2d ret("band_lim_wavelength",size(band_lims_wvn,1),size(band_lims_wvn,2));
    // for (int j = 1; j <= size(band_lims_wvn,2); j++) {
    //   for (int i = 1; i <= size(band_lims_wvn,1); i++) {
    parallel_for_cpu_serial( Bounds<2>( size(band_lims_wvn,2) , size(band_lims_wvn,1) ) , YAKL_LAMBDA (int j, int i) {
      if (this->is_initialized()) {
        ret(i,j) = 1._wp / this->band_lims_wvn(i,j);
      } else {
        ret(i,j) = 0._wp;
      }
    });
    return ret;
  }


  // Are the bands of two objects the same? (same number, same wavelength limits)
  bool bands_are_equal(OpticalProps const &rhs) {
    if ( this->get_nband() != rhs.get_nband() || this->get_nband() == 0) { return false; }
    // for (int j=1 ; j <= size(this->band_lims_wvn,2); j++) {
    //   for (int i=1 ; i <= size(this->band_lims_wvn,1); i++) {
    parallel_for_cpu_serial( Bounds<2>( size(this->band_lims_wvn,2) , size(this->band_lims_wvn,1) ) , YAKL_LAMBDA (int j, int i) {
      if ( abs( this->band_lims_wvn(i,j) - rhs.band_lims_wvn(i,j) ) > 5*epsilon(this->band_lims_wvn) ) {
        return false;
      }
    });
    return true;
  }


  // Is the g-point structure of two objects the same?
  //   (same bands, same number of g-points, same mapping between bands and g-points)
  bool gpoints_are_equal(OpticalProps const &rhs) {
    if ( ! this->bands_are_equal(rhs) || this->get_ngpt() != rhs.get_ngpt() ) { return false; }
    // for (int i=1; i <= size(this->gpt2bnd,1); i++) {
    parallel_for_cpu_serial( Bounds<1>(size(this->gpt2bnd,1)) , YAKL_LAMBDA (int i) {
      if ( this->gpt2band(i) != rhs.gpt2band(i) ) { return false; }
    });
    return true;
  }


  // Expand an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  real1d expand(real1d &arr_in) {
    real1d ret("arr_out",size(this->gpt2band,1));
    // do iband=1,this%get_nband()
    // TODO: I don't know if this needs to be serialize or not at first glance. Need to look at it more.
    parallel_for_cpu_serial( Bounds<1>(1) , YAKL_LAMBDA (int dummy) {
      for (int iband = 1 ; iband <= this->get_nband() ; iband++) {
        for (int i=this->band2gpt(1,iband) ; i <= this->band2gpt(2,iband) ; i++) {
          ret(i) = arr_in(iband);
        }
      }
    });
  }


  void set_name( std::string name ) { this->name = name; }


  void get_name() { return this->name; }
};



class OpticalPropsArry : public OpticalProps {
public:
  real3d tau; // optical depth (ncol, nlay, ngpt)

  YAKL_INLINE int get_ncol() { if (allocated(tau)) { return size(this->tau,1); } else { return 0; } }
  YAKL_INLINE int get_nlay() { if (allocated(tau)) { return size(this->tau,2); } else { return 0; } }
};



// Not implementing get_subset because it isn't used
class OpticalProps1scl : public OpticalPropsArry {
public:
  void validate() {
    if (! allocated(this->tau)) { stoprun("validate: tau not allocated/initialized"); }
    if (anyLT(this->tau,0._wp)) { stoprun("validate: tau values out of range"); }
  }


  void delta_scale(real3d const &dummy) { }


  void alloc_1scl(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("OpticalProps1scl::alloc_1scl: spectral discretization hasn't been provided"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("OpticalProps1scl::alloc_1scl: must provide > 0 extents for ncol, nlay"); }
    this->tau = real3d("tau",ncol,nlay,this->get_ngpt());
  }


  // Initialization by specifying band limits and possibly g-point/band mapping
  void alloc_1scl(int ncol, int nlay, real2d &band_lims_wvn, int2d const &band_lims_gpt=int2d(), std::string name="") {
    this->init(band_lims_wvn, band_lims_gpt, name);
    this->alloc_1scl(ncol, nlay);
  }


  void alloc_1scl(int ncol, int nlay, OpticalProps const &opIn, std::string name="") {
    if (this->is_initialized()) { this->finalize(); }
    this->init(opIn.get_band_lims_wavenumber(), opIn.get_band_lims_gpoint(), name);
    this->alloc_1scl(ncol, nlay);
  }


  void increment(OpticalProps1scl &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol()
    int nlay = that.get_nlay()
    int ngpt = that.get_ngpt()
    if (this->gpoints_are_equal(that)) {
      increment_1scalar_by_1scalar(ncol, nlay, ngpt, that.tau, this->tau);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, that.tau, this->tau, that.get_nband(), that.get_band_lims_gpoint());
    }
  }


  void increment(OpticalProps2str &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol()
    int nlay = that.get_nlay()
    int ngpt = that.get_ngpt()
    if (this->gpoints_are_equal(that)) {
      increment_1scalar_by_2stream(ncol, nlay, ngpt, that.tau, this->tau, this->ssa);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, that.tau, this->tau, this->ssa, that.get_nband(), that.get_band_lims_gpoint());
    }
  }
};



// Not implementing get_subset because it isn't used
class OpticalProps2str : public OpticalPropsArry {
public:
  real3d ssa // single-scattering albedo (ncol, nlay, ngpt)
  real3d g   // asymmetry parameter (ncol, nlay, ngpt)


  void validate() {
    if ( ! allocated(this%tau) || ! allocated(this%ssa) || ! allocated(this%g) ) {
      stoprun("validate: arrays not allocated/initialized");
    }
    int d1 = size(this->tau,1);
    int d2 = size(this->tau,2);
    int d3 = size(this->tau,3);
    if ( d1 != size(this->ssa,1) || d2 != size(this->ssa,2) || d3 != size(this->ssa,3) || 
         d1 != size(this->g  ,1) || d2 != size(this->g  ,2) || d3 != size(this->g  ,3) {
      stoprun("validate: arrays not sized consistently");
    }
    if (anyLT(this->tau,0._wp)                          ) { stoprun("validate: tau values out of range"); }
    if (anyLT(this->ssa,0._wp) || anyGT(this->ssa,1._wp)) { stoprun("validate: ssa values out of range"); }
    if (anyLT(this->g  ,0._wp) || anyGT(this->g  ,1._wp)) { stoprun("validate: g   values out of range"); }
  }


  void delta_scale(real3d const &forward=real3d()) {
    // Forward scattering fraction; g**2 if not provided
    int ncol = this->get_ncol();
    int nlay = this->get_nlay();
    int ngpt = this->get_ngpt();
    if (allocated(forward)) {
      if (size(forward,1) != ncol || size(forward,2) != nlay || size(forward,3) != ngpt) {
        stoprun("delta_scale: dimension of 'forward' don't match optical properties arrays");
      }
      if (anyLT(forward,0._wp) || anyGT(forward,1._wp)) { stoprun("delta_scale: values of 'forward' out of bounds [0,1]"); }
      delta_scale_2str_kernel(ncol, nlay, ngpt, this->tau, this->ssa, this->g, forward)
    } else {
      delta_scale_2str_kernel(ncol, nlay, ngpt, this->tau, this->ssa, this->g)
    }
  }


  void alloc_2str(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("optical_props%alloc: spectral discretization hasn't been provided"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("optical_props%alloc: must provide positive extents for ncol, nlay"); }
    this->tau = real3d("tau",ncol,nlay,this->get_ngpt());
    this->ssa = real3d("ssa",ncol,nlay,this->get_ngpt());
    this->g   = real3d("g  ",ncol,nlay,this->get_ngpt());
  }


  void alloc_2str(int ncol, int nlay, real2d const &band_lims_wvn, int2d const &band_lims_gpt=int2d(), std::string name="") {
    this->init(band_lims_wvn, band_lims_gpt, name);
    this->alloc_2str(ncol, nlay);
  }


  void alloc_2str(int ncol, int nlay, OpticalProps const &opIn, std::string name="") {
    if (this->is_initialized()) { this->finalize(); }
    this->init(opIn.get_band_lims_wavenumber(), opIn.get_band_lims_gpoint(), name);
    this->alloc_2str(ncol, nlay);
  }


  void increment(OpticalProps1scl &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol()
    int nlay = that.get_nlay()
    int ngpt = that.get_ngpt()
    if (this->gpoints_are_equal(that)) {
      increment_2stream_by_1scalar(ncol, nlay, ngpt, that.tau, that.ssa, this->tau);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, that.tau, that.ssa, this->tau, that.get_nband(), that.get_band_lims_gpoint());
    }
  }


  void increment(OpticalProps2str &that) {
    if (! this->bands_are_equal(that)) { stoprun("OpticalProps::increment: optical properties objects have different band structures"); }
    int ncol = that.get_ncol()
    int nlay = that.get_nlay()
    int ngpt = that.get_ngpt()
    if (this->gpoints_are_equal(that)) {
      increment_2stream_by_2stream(ncol, nlay, ngpt, that.tau, that.ssa, that.g, this->tau, this->ssa, this->g);
    } else {
      if (this->get_ngpt() != that.get_nband()) {
        stoprun("OpticalProps::increment: optical properties objects have incompatible g-point structures");
      }
      inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt, that.tau, that.ssa, that.g, this->tau, this->ssa, this->g, that.get_nband(), that.get_band_lims_gpoint());
    }
  }
};



