
#include "mo_rrtmgp_util_reorder_kernels.h"

extern "C" void reorder_123x321_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p) {
  umgReal3d array_in ( "array_in"  , array_in_p  , d3 , d2 , d1 );
  umgReal3d array_out( "array_out" , array_out_p , d1 , d2 , d3 );

  for (int i3=0; i3<d3; i3++) {
    for (int i2=0; i2<d2; i2++) {
      for (int i1=0; i1<d1; i1++) {
        array_out(i1,i2,i3) = array_in(i3,i2,i1);
      }
    }
  }
}

