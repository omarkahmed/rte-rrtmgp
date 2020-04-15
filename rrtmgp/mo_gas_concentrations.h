
#pragma once

#include "const.h"

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
  real3d   concs;
  int      ncol;
  int      nlay;
  int      ngas;


  GasConcs() {
    ncol = 0;
    nlay = 0;
  }


  ~GasConcs() {
    reset();
  }


  void reset() {
    gas_name = string1d();
    concs    = real3d();
    ncol = 0;
    nlay = 0;
  }


  void init(string1d &gas_names , int ncol , int nlay) {
    this->ngas = size(gas_names,1);

    // Transform gas names to lower case, and check for empty strings
    for (int i=1; i<=ngas; i++) {
      gas_names(i) = lower_case( gas_names(i) );
      // Empty string
      if (gas_names(i) == "") { stoprun("ERROR: GasConcs::init(): must provide non-empty gas names"); }
      // Duplicate gas name
      for (int j=i+1; j<=ngas; j++) {
        if ( gas_names(i) == gas_names(j) ) { stoprun("GasConcs::init(): duplicate gas names aren't allowed"); }
      }
    }
    
    // Allocate fixed-size arrays
    this->reset();
    this->ncol = ncol;
    this->nlay = nlay;
    this->gas_name = string1d("gas_name",ngas);
    this->concs    = real3d  ("concs"   ,ncol,nlay,ngas);

    for (int i=1; i<=ngas; i++) {
      this->gas_name(i) = gas_names(i);
    }
  }


  // Set concentrations --- scalar, 1D, 2D
  void set_vmr(std::string gas, real w) {
    if (w < 0._wp || w > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name was not provided at initialization");
    }
    for (int ilay=1; ilay<=this->nlay; ilay++) {
      for (int icol=1; icol<=this->ncol; icol++) {
        this->concs(icol,ilay,igas) = w;
      }
    }
  }


  void set_vmr(std::string gas, real1d &w) {
    if (size(w,1) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)"); }
    for (int i=1; i<=size(w,1); i++) {
      if (w(i) < 0._wp || w(i) > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
    }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization");
    }
    for (int ilay=1; ilay<=this->nlay; ilay++) {
      for (int icol=1; icol<=this->ncol; icol++) {
        this->concs(icol,ilay,igas) = w(ilay);
      }
    }
  }
  

  void set_vmr(std::string gas, real2d &w) {
    if (size(w,1) != this->ncol) { stoprun("GasConcs::set_vmr: different dimension (ncol)" ); }
    if (size(w,2) != this->nlay) { stoprun("GasConcs::set_vmr: different dimension (nlay)" ); }
    for (int j=1; j<=size(w,2); j++) {
      for (int i=1; i<=size(w,1); i++) {
        if (w(i,j) < 0._wp || w(i,j) > 1._wp) { stoprun("GasConcs::set_vmr(): concentrations should be >= 0, <= 1"); }
      }
    }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) {
      stoprun("GasConcs::set_vmr(): trying to set a gas whose name not provided at initialization" );
    }
    for (int ilay=1; ilay<=this->nlay; ilay++) {
      for (int icol=1; icol<=this->ncol; icol++) {
        this->concs(icol,ilay,igas) = w(icol,ilay);
      }
    }
  }


  // 2D array (col, lay)
  void get_vmr(std::string gas, real2d &array) {
    if (this->ncol != size(array,1)) { stoprun("ty_gas_concs->get_vmr; gas array is wrong size (ncol)" ); }
    if (this->nlay != size(array,2)) { stoprun("ty_gas_concs->get_vmr; gas array is wrong size (nlay)" ); }
    int igas = this->find_gas(gas);
    if (igas == GAS_NOT_IN_LIST) { stoprun("GasConcs::get_vmr; gas not found" ); }
    for (int ilay=1; ilay<=size(array,2); ilay++) {
      for (int icol=1; icol<=size(array,1); icol++) {
        array(icol,ilay) = this->concs(icol,ilay,igas);
      }
    }
  }


  int get_num_gases() { return size(gas_name,1); }
  

  string1d get_gas_names() { return gas_name; }


  // find gas in list; GAS_NOT_IN_LIST if not found
  int find_gas(std::string gas) {
    if (ngas == 0) { return 0; }
    for (int igas=1; igas<=ngas; igas++) {
      if ( lower_case(gas) == this->gas_name(igas) ) {
        return igas;
      }
    }
    return GAS_NOT_IN_LIST;
  }


};


