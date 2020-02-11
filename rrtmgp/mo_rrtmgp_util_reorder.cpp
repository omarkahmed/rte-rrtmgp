
#include "mo_rrtmgp_util_reorder.h"
#include "mo_rrtmgp_util_reorder_kernels.h"

// (x,y,z) -> (z,x,y)
extern "C" void reorder123x312(int d1, int d2, int d3, real *array_p, real *array_out_p) {
  umgReal3d array    ("array    ",array_p    ,d1,d2,d3);
  umgReal3d array_out("array_out",array_out_p,d3,d1,d2);

  reorder_123x312_kernel(d1, d2, d3, array, array_out);
  std::cout << "WARNING: THIS IS NOT TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// (x,y,z) -> (z,y,x)
extern "C" void reorder123x321(int d1, int d2, int d3, real *array_p, real *array_out_p) {
  umgReal3d array    ("array    ",array_p    ,d1,d2,d3);
  umgReal3d array_out("array_out",array_out_p,d3,d2,d1);

  reorder_123x321_kernel(d1, d2, d3, array, array_out);
}



