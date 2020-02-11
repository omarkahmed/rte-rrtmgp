
#pragma once

#include "const.h"

extern "C" void reorder_123x321_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p);

extern "C" void reorder_123x312_kernel(int d1, int d2, int d3, real *array_in_p, real *array_out_p);

