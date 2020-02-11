
#include "mo_rrtmgp_util_reorder_kernels.h"

extern "C" void reorder_123x321_kernel(int d1, int d2, int d3, real * __restrict__ array_in_p, real * __restrict__ array_out_p) {
  umgReal3d array_in ( "array_in"  , array_in_p  , d1 , d2 , d3 );
  umgReal3d array_out( "array_out" , array_out_p , d3 , d2 , d1 );

  parallel_for_cpu_serial( Bounds<3>({1,d3},{1,d2},{1,d1}) , YAKL_LAMBDA ( int const indices[] ) {
    int i3, i2, i1;
    yakl::storeIndices( indices , i3,i2,i1 );

    array_out(i3,i2,i1) = array_in(i1,i2,i3);
  });
}

