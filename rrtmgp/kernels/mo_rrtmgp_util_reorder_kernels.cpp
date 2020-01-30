
#include "mo_rrtmgp_util_reorder_kernels.h"

extern "C" void reorder_123x321_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p) {
  umgReal3d array_in ( "array_in"  , array_in_p  , d3 , d2 , d1 );
  umgReal3d array_out( "array_out" , array_out_p , d1 , d2 , d3 );

  // for (int k=0; k<d3; k++) {
  //   for (int j=0; j<d2; j++) {
  //     for (int i=0; i<d1; i++) {
  yakl::parallel_for( d3*d2*d1 , YAKL_LAMBDA (int iGlob) {
    int i3, i2, i1;
    yakl::unpackIndices( iGlob , d3,d2,d1 , i3,i2,i1 );

    array_out(i1,i2,i3) = array_in(i3,i2,i1);
  });
}

