
#pragma once

#include "const.h"

using yakl::memHost;
using yakl::memDevice;

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
// Encapsulates a collection of volume mixing ratios (concentrations) of gases.
//   Each concentration is associated with a name, normally the chemical formula.
//
// Values may be provided as scalars, 1-dimensional profiles (nlay), or 2-D fields (ncol,nlay).
//   (nlay and ncol are determined from the input arrays; self-consistency is enforced)
//   example:
//   error_msg = gas_concs.set_vmr('h2o', values(:,:))
//   error_msg = gas_concs.set_vmr('o3' , values(:)  )
//   error_msg = gas_concs.set_vmr('co2', value      )
//
// Values can be requested as profiles (valid only if there are no 2D fields present in the object)
//   or as 2D fields. Values for all columns are returned although the entire collection
//   can be subsetted in the column dimension
//
// Subsets can be extracted in the column dimension
//
// Functions return strings. Non-empty strings indicate an error.
//
// -------------------------------------------------------------------------------------------------


class GasConcs {
public:
  static int constexpr GAS_NOT_IN_LIST = -1;
  string1d gas_name;
  FArray<real2d,1,memDevice> concs;
  int ncol;
  int nlay;


  GasConcs() {
    ncol = 0;
    nlay = 0;
  }


  ~GasConcs() {
    reset();
  }


  void reset() {
    gas_name = FArray<std::string,1,memHost  >();
    concs    = FArray<real2d     ,1,memDevice>();
    ncol = 0;
    nlay = 0;
  }


  void init(FArray<std::string,1,memHost> &gas_names) {
    int ngas = size(gas_names,1);

    // Transform gas names to lower case, and check for empty strings
    for (int i=1; i<=ngas; i++) {
      gas_names(i) = lower_case( gas_names(i) );
      if (gas_names(i).length() == 0) { stoprun("ERROR: GasConcs::init(): must provide non-empty gas names"); }
    }

    // Check for duplicate gase names
    for (int i=1; i<=ngas-1; i++) {
      for (int j=i+1; j<=ngas; j++) {
        if ( gas_names(i) == gas_names(j) ) { stoprun("GasConcs::init(): duplicate gas names aren't allowed"); }
      }
    }
    
    // Allocate fixed-size arrays
    this->reset();
    this->gas_name = FArray<std::string,1,memHost>  ("gas_name",ngas);
    this->concs    = FArray<real2d     ,1,memDevice>("concs"   ,ngas);

    for (int i=1; i<=ngas; i++) {
      this->gas_name(i) = gas_names(i);
    }
  }


  // Set concentrations --- scalar, 1D, 2D
  void set_vmr(std::string gas, real w) {
    int igas = this->find_gas(gas);

    if (w < 0._wp || w > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name was not provided at initialization");
    }
    
    this->concs(igas) = real2d("concs",1,1);
    this->concs(igas)(1,1) = w;
  }


  void set_vmr(std::string gas, real1d &w) {
    int igas = this->find_gas(gas);

    for (int i=1; i<=size(w,1); i++) {
      if (w(i) < 0._wp || w(i) > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    }
    if (this->nlay > 0 && size(w,1) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)"); }
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization");
    }

    this->nlay = size(w,1);
    this->concs(igas) = real2d("concs",1,this->nlay);
    for (int ilay=1; ilay<=this->nlay; ilay++) {
      this->concs(igas)(1,ilay) = w(ilay);
    }
  }
  

  void set_vmr(std::string gas, real2d &w) {
    int igas = this->find_gas(gas);

    for (int j=1; j<=size(w,2); j++) {
      for (int i=1; i<=size(w,1); i++) {
        if (w(i,j) < 0._wp || w(i,j) > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
      }
    }
    if (this->ncol > 0 && size(w,1) != this->ncol) { stoprun("GasConcs::set_vmr: different dimension (ncol)" ); }
    if (this->nlay > 0 && size(w,2) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)" ); }
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization" );
    }
    
    this->ncol = size(w,1);
    this->nlay = size(w,2);
    this->concs(igas) = real2d("concs",this->ncol,this->nlay);

    for (int ilay=1; ilay<=this->nlay; ilay++) {
      for (int icol=1; icol<=this->ncol; icol++) {
        this->concs(igas)(icol,ilay) = w(icol,ilay);
      }
    }
  }


  // Return volume mixing ratio as 1D or 2D array
  // 1D array ( lay depdendence only)
  void get_vmr(std::string gas, real1d &array) {
    int igas = this->find_gas(gas);

    if (igas == GAS_NOT_IN_LIST) { stoprun("GasConcs::get_vmr; gas not found" ); }
    if (this->concs(igas).totElems() == 0 ) {
      stoprun("GasConcs::get_vmr; gas concentration hasn't been set" );
    }
    if (size(this->concs(igas),1) > 1) { // Are we requesting a single profile when many are present?
      stoprun("GasConcs::get_vmr; gas requesting single profile but many are available" );
    }
    if (this->nlay > 0 && this->nlay != size(array,1)) {
      stoprun("GasConcs::get_vmr; gas array is wrong size (nlay)" );
    }

    if ( size(this->concs(igas),2) > 1) {
      for (int i=1; i<=size(array,1); i++) {
        array(i) = this->concs(igas)(1,i);
      }
    } else {
      for (int i=1; i<=size(array,1); i++) {
        array(i) = this->concs(igas)(1,1);
      }
    }
  }


  // 2D array (col, lay)
  void get_vmr(std::string gas, real2d &array) {
    int igas = this->find_gas(gas);

    if (igas == GAS_NOT_IN_LIST) { stoprun("GasConcs::get_vmr; gas not found" ); }
    if (this->concs(igas).totElems() == 0) {
      stoprun("GasConcs::get_vmr; gas concentration hasn't been set" );
    }
    if (this->ncol > 0 && this->ncol != size(array,1)) {
      stoprun("ty_gas_concs->get_vmr; gas array is wrong size (ncol)" );
    }
    if (this->nlay > 0 && this->nlay != size(array,2)) {
      stoprun("ty_gas_concs->get_vmr; gas array is wrong size (nlay)" );
    }

    if ( size(this->concs(igas),1) > 1) {       // Concentration stored as 2D
      for (int ilay=1; ilay<=size(array,2); ilay++) {
        for (int icol=1; icol<=size(array,1); icol++) {
          array(icol,ilay) = this->concs(igas)(icol,ilay);
        }
      }
    } else if ( size(this->concs(igas),2) > 1) { // Concentration stored as 1D
      for (int ilay=1; ilay<=size(array,2); ilay++) {
        for (int icol=1; icol<=size(array,1); icol++) {
         array(icol,ilay) = this->concs(igas)(1,ilay);
        }
      }
    } else {                                         // Concentration stored as scalar
      for (int ilay=1; ilay<=size(array,2); ilay++) {
        for (int icol=1; icol<=size(array,1); icol++) {
          array(icol,ilay) = this->concs(igas)(1,1);
        }
      }
    }
  }


  // Extract a subset of n columns starting with column 'start'
  void get_subset_range(int start, int n, GasConcs &subset) {
    if ( n <= 0    ) { stoprun("GasConcs::get_vmr: Asking for 0 or fewer columns "   ); }
    if ( start < 1 ) { stoprun("GasConcs::get_vmr: Asking for columns outside range" ); }
    if ( this->ncol > 0 && start > this->ncol || start+n-1 > this->ncol ) {
      stoprun("GasConcs::get_vmr: Asking for columns outside range" );
    }

    subset.reset();

    int ngas = size(this->gas_name,1);
    subset.gas_name = FArray<std::string,1,memHost  >("gas_name",ngas);
    subset.concs    = FArray<real2d     ,1,memDevice>("concs"   ,ngas);
    subset.nlay = this->nlay;
    subset.ncol = this->ncol > 0 ? n : 0;
    for (int i=1; i<=ngas; i++) {
      subset.gas_name(i) = this->gas_name(i);

      // Preserve scalar/1D/2D representation in subset,
      //   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
      int ncol_subset = min( max( subset.ncol , 1 ) , size(this->concs(i),1) );
      int nlay_subset = min(      subset.nlay       , size(this->concs(i),2) );
      subset.concs(i) = real2d("concs",ncol_subset,nlay_subset);
      if (size(this->concs(i),1) > 1) {      // Concentration stored as 2D
        for (int ilay=1; ilay<=nlay_subset; ilay++) {
          for (int icol=1; icol<=ncol_subset; icol++) {
            subset.concs(i)(icol,ilay) = this->concs(i)(icol+start-1,ilay);
          }
        }
      } else {
        for (int ilay=1; ilay<=size(this->concs(i),2); ilay++) {
          for (int icol=1; icol<=size(this->concs(i),1); icol++) {
            subset.concs(i)(icol,ilay) = this->concs(i)(icol,ilay);
          }
        }
      }
    }
  }


  int get_num_gases() { return size(gas_name,1); }
  

  FArray<std::string,1,memHost> get_gas_names() { return gas_name; }


  // find gas in list; GAS_NOT_IN_LIST if not found
  int find_gas(std::string gas) {
    int ngas = size(this->gas_name,1);
    if (ngas == 0) { return 0; }
    for (int igas=1; igas<=ngas; igas++) {
      if ( lower_case(gas) == this->gas_name(igas) ) {
        return igas;
      }
    }
    return GAS_NOT_IN_LIST;
  }


};


