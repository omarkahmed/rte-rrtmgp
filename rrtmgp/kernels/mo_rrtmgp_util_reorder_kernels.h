
#pragma once

#include "const.h"


YAKL_INLINE void reorder_123x321_kernel(int d1, int d2, int d3, real3d const &array_in, real3d &array_out) {
  parallel_for_cpu_serial( Bounds<3>({1,d3},{1,d2},{1,d1}) , YAKL_LAMBDA ( int const indices[] ) {
    int i3, i2, i1;
    yakl::storeIndices( indices , i3,i2,i1 );

    array_out(i3,i2,i1) = array_in(i1,i2,i3);
  });
}


YAKL_INLINE void reorder_123x312_kernel(int d1, int d2, int d3, real3d const &array_in, real3d &array_out) {
  // do i3 = 1 , d3
  //   do i2 = 1 , d2
  //     do i1 = 1 , d1
  parallel_for( Bounds<3>({1,d3},{1,d2},{1,d1}) , YAKL_LAMBDA (int const indices[]) {
    int i3, i2, i1;
    storeIndices( indices , i3, i2, i1 );
    array_out(i3,i1,i2) = array_in(i1,i2,i3);
  });
}


