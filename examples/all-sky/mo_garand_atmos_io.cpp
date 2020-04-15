
#include "mo_garand_atmos_io.h"

// Read data from the input file
void read_atmos(std::string input_file, realHost2d &p_lay, realHost2d &t_lay, realHost2d &p_lev, realHost2d &t_lev,
                GasConcs &gas_concs, realHost2d &col_dry) {
  yakl::SimpleNetCDF io;
  io.open(input_file,yakl::MODE_READ);

  int ncol = io.getDimSize("col");
  int nlay = io.getDimSize("lay");

  io.read(p_lay,"p_lay");
  io.read(t_lay,"t_lay");
  io.read(p_lev,"p_lev");
  io.read(t_lev,"t_lev");

  int ngas = 8;
  string1d gas_names("gas_names",ngas);
  gas_names(1) = std::string("h2o");
  gas_names(2) = std::string("co2");
  gas_names(3) = std::string("o3" );
  gas_names(4) = std::string("n2o");
  gas_names(5) = std::string("co" );
  gas_names(6) = std::string("ch4");
  gas_names(7) = std::string("o2" );
  gas_names(8) = std::string("n2" );

  gas_concs.init(gas_names,ncol,nlay);

  for (int igas=1 ; igas <= ngas ; igas++) {
    std::string vmr_name = "vmr_"+gas_names(igas);
    if ( ! io.varExists(vmr_name) ) { throw "ERROR: gas does not exist in input flie"; }
    real2d tmp;
    io.read(tmp,vmr_name);  // Automatically copies it to the device
    gas_concs.set_vmr( gas_names(igas) , tmp );
  }

  if ( io.varExists("col_dry") ) { io.read(col_dry,"col_dry"); }

  io.close();
}


