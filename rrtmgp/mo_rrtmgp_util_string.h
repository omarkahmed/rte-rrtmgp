
#pragma once


inline std::string lower_case( std::string in ) {
  std::for_each( in.begin() , in.end() , [] (char & c) { c = ::tolower(c); } );
  return in;
}


inline bool string_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return true; }
  }
  return false;
}


inline int string_loc_in_array(std::string str, string1d const &arr) {
  for (int i=1; i <= size(arr); i++) {
    if ( lower_case(str) == lower_case(arr(i)) ) { return i; }
  }
  return -1;
}


