
#pragma once

#include "YAKL.h"
#include "FArray.h"
#include "SArray.h"
#include "FSArray.h"
#include <iostream>

typedef double real;

YAKL_INLINE real constexpr operator"" _wp( long double x ) {
  return static_cast<real>(x);
}


using yakl::fortran::parallel_for;
using yakl::fortran::Bounds;
using yakl::storeIndices;


typedef yakl::FArray<real,1,yakl::memDevice> umgReal1d;
typedef yakl::FArray<real,2,yakl::memDevice> umgReal2d;
typedef yakl::FArray<real,3,yakl::memDevice> umgReal3d;
typedef yakl::FArray<real,4,yakl::memDevice> umgReal4d;
typedef yakl::FArray<real,5,yakl::memDevice> umgReal5d;
typedef yakl::FArray<real,6,yakl::memDevice> umgReal6d;
typedef yakl::FArray<real,7,yakl::memDevice> umgReal7d;

typedef yakl::FArray<real,1,yakl::memDevice> real1d;
typedef yakl::FArray<real,2,yakl::memDevice> real2d;
typedef yakl::FArray<real,3,yakl::memDevice> real3d;
typedef yakl::FArray<real,4,yakl::memDevice> real4d;
typedef yakl::FArray<real,5,yakl::memDevice> real5d;
typedef yakl::FArray<real,6,yakl::memDevice> real6d;
typedef yakl::FArray<real,7,yakl::memDevice> real7d;

typedef yakl::FArray<real,1,yakl::memHost> realHost1d;
typedef yakl::FArray<real,2,yakl::memHost> realHost2d;
typedef yakl::FArray<real,3,yakl::memHost> realHost3d;
typedef yakl::FArray<real,4,yakl::memHost> realHost4d;
typedef yakl::FArray<real,5,yakl::memHost> realHost5d;
typedef yakl::FArray<real,6,yakl::memHost> realHost6d;
typedef yakl::FArray<real,7,yakl::memHost> realHost7d;

typedef yakl::FArray<int,1,yakl::memDevice> umgInt1d;
typedef yakl::FArray<int,2,yakl::memDevice> umgInt2d;
typedef yakl::FArray<int,3,yakl::memDevice> umgInt3d;
typedef yakl::FArray<int,4,yakl::memDevice> umgInt4d;
typedef yakl::FArray<int,5,yakl::memDevice> umgInt5d;
typedef yakl::FArray<int,6,yakl::memDevice> umgInt6d;
typedef yakl::FArray<int,7,yakl::memDevice> umgInt7d;

typedef yakl::FArray<int,1,yakl::memDevice> int1d;
typedef yakl::FArray<int,2,yakl::memDevice> int2d;
typedef yakl::FArray<int,3,yakl::memDevice> int3d;
typedef yakl::FArray<int,4,yakl::memDevice> int4d;
typedef yakl::FArray<int,5,yakl::memDevice> int5d;
typedef yakl::FArray<int,6,yakl::memDevice> int6d;
typedef yakl::FArray<int,7,yakl::memDevice> int7d;

typedef yakl::FArray<int,1,yakl::memHost> intHost1d;
typedef yakl::FArray<int,2,yakl::memHost> intHost2d;
typedef yakl::FArray<int,3,yakl::memHost> intHost3d;
typedef yakl::FArray<int,4,yakl::memHost> intHost4d;
typedef yakl::FArray<int,5,yakl::memHost> intHost5d;
typedef yakl::FArray<int,6,yakl::memHost> intHost6d;
typedef yakl::FArray<int,7,yakl::memHost> intHost7d;

typedef yakl::FArray<bool,1,yakl::memDevice> umgBool1d;
typedef yakl::FArray<bool,2,yakl::memDevice> umgBool2d;
typedef yakl::FArray<bool,3,yakl::memDevice> umgBool3d;
typedef yakl::FArray<bool,4,yakl::memDevice> umgBool4d;
typedef yakl::FArray<bool,5,yakl::memDevice> umgBool5d;
typedef yakl::FArray<bool,6,yakl::memDevice> umgBool6d;
typedef yakl::FArray<bool,7,yakl::memDevice> umgBool7d;

typedef yakl::FArray<bool,1,yakl::memDevice> bool1d;
typedef yakl::FArray<bool,2,yakl::memDevice> bool2d;
typedef yakl::FArray<bool,3,yakl::memDevice> bool3d;
typedef yakl::FArray<bool,4,yakl::memDevice> bool4d;
typedef yakl::FArray<bool,5,yakl::memDevice> bool5d;
typedef yakl::FArray<bool,6,yakl::memDevice> bool6d;
typedef yakl::FArray<bool,7,yakl::memDevice> bool7d;

typedef yakl::FArray<bool,1,yakl::memHost> boolHost1d;
typedef yakl::FArray<bool,2,yakl::memHost> boolHost2d;
typedef yakl::FArray<bool,3,yakl::memHost> boolHost3d;
typedef yakl::FArray<bool,4,yakl::memHost> boolHost4d;
typedef yakl::FArray<bool,5,yakl::memHost> boolHost5d;
typedef yakl::FArray<bool,6,yakl::memHost> boolHost6d;
typedef yakl::FArray<bool,7,yakl::memHost> boolHost7d;


template <class T> YAKL_INLINE  T max(T v1, T v2) { return v1 > v2 ? v1 : v2; }
template <class T> YAKL_INLINE  T min(T v1, T v2) { return v1 < v2 ? v1 : v2; }

template <class T> YAKL_INLINE  T merge(T t, T f, bool cond) { return cond ? t : f; }

template <class T> YAKL_INLINE  int size(T arr, int dim) { return arr.dimension[dim-1]; }


inline void stoprun( std::string str ) {
  throw str;
}

inline std::string lower_case( std::string str ) {
  std::string ret = str;
  std::transform(ret.begin(), ret.end(), ret.begin(), [](unsigned char c){ return std::tolower(c); });
  return ret;
}



