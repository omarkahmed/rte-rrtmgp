#!/usr/bin/env bash

source ${HOME}/e3sm-comps/rte-rrtmgp/cpp/test/intel.sh

module load intel-comp-rt/agama-ci-prerelease/398

./cmakeclean.sh
./cmakescript.sh
cd build
make -j 32

export LD_LIBRARY_PATH=$HOME/e3sm-libs/oneapi-ifort/packages/netcdf-serial/lib:$HOME/e3sm-libs/oneapi-ifort/packages/hdf5-serial/lib:$HOME/e3sm-libs/oneapi-ifort/packages/szip/lib:$LD_LIBRARY_PATH

ctest
