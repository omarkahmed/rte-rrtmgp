
#include "mo_rrtmgp_util_reorder_kernels.h"



extern "C" void reorder_123x321_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p) {
  umgReal3d array_in ( "array_in"  , array_in_p  , d1 , d2 , d3 );
  umgReal3d array_out( "array_out" , array_out_p , d3 , d2 , d1 );

  parallel_for_cpu_serial( Bounds<3>({1,d3},{1,d2},{1,d1}) , YAKL_LAMBDA ( int const indices[] ) {
    int i3, i2, i1;
    yakl::storeIndices( indices , i3,i2,i1 );

    array_out(i3,i2,i1) = array_in(i1,i2,i3);
  });
}



extern "C" void reorder_123x312_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p) {
  umgReal3d array_in ( "array_in"  , array_in_p  , d1 , d2 , d3 );
  umgReal3d array_out( "array_out" , array_out_p , d3 , d1 , d2 );

  // do i3 = 1 , d3
  //   do i2 = 1 , d2
  //     do i1 = 1 , d1
  parallel_for( Bounds<3>({1,d3},{1,d2},{1,d1}) , YAKL_LAMBDA (int const indices[]) {
    int i3, i2, i1;
    storeIndices( indices , i3, i2, i1 );
    array_out(i3,i1,i2) = array_in(i1,i2,i3);
  });
  std::cout << "WARNING: THIS IS NOT TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}

