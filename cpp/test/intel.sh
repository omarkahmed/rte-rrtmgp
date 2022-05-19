#!/bin/bash
module purge
module load intel-nightly/20220516
module load intel/mpi/2021.5.0

export ARCH="SYCL"
export CXXFLAGS="-fsycl -sycl-target=spirv64_gen"
export NCINCLUDE="`$HOME/e3sm-libs/oneapi-ifort/packages/netcdf-serial/bin/nc-config --includedir`"
export NCFLAGS="`$HOME/e3sm-libs/oneapi-ifort/packages/netcdf-serial/bin/nc-config --libs`"
export CC=icx
export CXX=icpx
#export CXX=dpcpp
export FC=ifx
export YAKLHOME="$HOME/e3sm-comps/YAKL"
