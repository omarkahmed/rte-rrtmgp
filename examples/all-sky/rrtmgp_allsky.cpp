
#include <iostream>
#include <cstdlib>
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "mo_gas_optics_rrtmgp.h"
#include "const.h"

int main(int argc , char **argv) {
  yakl::init();

  if (argc < 5) { stoprun("Error: Fewer than 4 command line arguments provided"); }
  std::string input_file        =      argv[1];
  std::string k_dist_file       =      argv[2];
  std::string cloud_optics_file =      argv[3];
  int ncol                      = atoi(argv[4]);
  int nloops = 1;
  if (argc >= 6) { nloops       = atoi(argv[5]); }
  if (ncol   <= 0) { stoprun("Error: Number of columns must be > 0"); }
  if (nloops <= 0) { stoprun("Error: Number of loops must be > 0"); }
  if (argc > 6) { std::cout << "WARNING: Using only 5 parameters. Ignoring the rest\n"; }
  if (input_file == "-h" || input_file == "--help") {
    std::cout << "./rrtmgp_allsky  input_file  absorption_coefficients_file  cloud_optics_file  ncol  [nloops]\n\n";
    exit(0);
  }

  // Read temperature, pressure, gas concentrations. Arrays are allocated as they are read
  realHost2d p_lay;
  realHost2d t_lay;
  realHost2d p_lev;
  realHost2d t_lev;
  GasConcs gas_concs;
  realHost2d col_dry;

  // Read data from the input file
  read_atmos(input_file, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

  int nlay = size(p_lay,2);

  // load data into classes
  GasOpticsRRTMGP k_dist;
  load_and_init(k_dist, k_dist_file, gas_concs);

  bool is_sw = k_dist.source_is_external();
  bool is_lw = ! is_sw;

  yakl::finalize();
}
