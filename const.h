
#pragma once

#include "YAKL.h"
#include "Array.h"
#include "SArray.h"

typedef double real;


typedef yakl::Array<real,yakl::memDevice> umgReal1d;
typedef yakl::Array<real,yakl::memDevice> umgReal2d;
typedef yakl::Array<real,yakl::memDevice> umgReal3d;
typedef yakl::Array<real,yakl::memDevice> umgReal4d;
typedef yakl::Array<real,yakl::memDevice> umgReal5d;
typedef yakl::Array<real,yakl::memDevice> umgReal6d;
typedef yakl::Array<real,yakl::memDevice> umgReal7d;

typedef yakl::Array<real,yakl::memDevice> real1d;
typedef yakl::Array<real,yakl::memDevice> real2d;
typedef yakl::Array<real,yakl::memDevice> real3d;
typedef yakl::Array<real,yakl::memDevice> real4d;
typedef yakl::Array<real,yakl::memDevice> real5d;
typedef yakl::Array<real,yakl::memDevice> real6d;
typedef yakl::Array<real,yakl::memDevice> real7d;

typedef yakl::Array<real,yakl::memHost> realHost1d;
typedef yakl::Array<real,yakl::memHost> realHost2d;
typedef yakl::Array<real,yakl::memHost> realHost3d;
typedef yakl::Array<real,yakl::memHost> realHost4d;
typedef yakl::Array<real,yakl::memHost> realHost5d;
typedef yakl::Array<real,yakl::memHost> realHost6d;
typedef yakl::Array<real,yakl::memHost> realHost7d;


template <class T> YAKL_INLINE  T max(T v1, T v2) { return v1 > v2 ? v1 : v2; }
template <class T> YAKL_INLINE  T min(T v1, T v2) { return v1 < v2 ? v1 : v2; }



