
#include "const.h"
#include <ctype.h>
#include <stdio.h>


extern "C" void lower_case( char *input_string , char *output_string ) {
  int i = 0;
  while (input_string[i] != '\0') {
    output_string[i] = tolower( input_string[i] );
    i++;
  }
  output_string[i] = '\0';
}
