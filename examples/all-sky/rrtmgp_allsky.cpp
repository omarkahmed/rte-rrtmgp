


int main(int argc , char **argv) {
  if (argc < 5) { throw "Error: Fewer than 4 command line arguments provided"; }
  std::string input_file        =      argv[1];
  std::string k_dist_file       =      argv[2];
  std::string cloud_optics_file =      argv[3];
  int ncol                      = atoi(argv[4]);
  int nloops = 1;
  if (argc >= 5) { nloops       = atoi(argv[5]); }
  if (ncol   <= 0) { throw "Error: Number of columns must be > 0"; }
  if (nloops <= 0) { throw "Error: Number of loops must be > 0"; }
  if (argc > 5) { std::cout << "WARNING: Using only 5 parameters. Ignoring the rest\n"; }

  // Read temperature, pressure, gas concentrations.
  //   Arrays are allocated as they are read
  call read_atmos(input_file, p_lay, t_lay, p_lev, t_lev, gas_concs_garand, col_dry)
}
