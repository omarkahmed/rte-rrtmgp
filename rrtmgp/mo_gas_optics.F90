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
! -------------------------------------------------------------------------------------------------
module mo_gas_optics
  use mo_rte_kind,           only: wp
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_arry

  type, extends(ty_optical_props), public :: ty_gas_optics
  end type
end module mo_gas_optics
