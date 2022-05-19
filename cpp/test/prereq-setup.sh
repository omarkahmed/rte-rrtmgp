#!/bin/bash

module load intel-nightly/20220516
module load intel/mpi/2021.5.0

export installation_home=$HOME/e3sm-libs/oneapi-ifort
export procs=32
export c_compiler=icc
export cxx_compiler=icpc
export fortran_compiler=ifort
export mpic_compiler=mpiicc
export mpicxx_compiler=mpiicpc
export mpifortran_compiler=mpiifort

mkdir $installation_home/sources
mkdir $installation_home/packages
cd $installation_home/sources
wget --inet4-only https://github.com/Kitware/CMake/releases/download/v3.22.2/cmake-3.22.2.tar.gz
#wget --inet4-only http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz  --no-check-certificate
wget --inet4-only https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
wget --inet4-only https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_10_6/source/hdf5-1.10.6.tar.gz
wget --inet4-only https://github.com/Unidata/netcdf-c/archive/v4.7.3.tar.gz
wget --inet4-only https://github.com/Unidata/netcdf-cxx4/archive/v4.3.1.tar.gz
wget --inet4-only https://github.com/Unidata/netcdf-fortran/archive/v4.4.4.tar.gz
wget --inet4-only https://parallel-netcdf.github.io/Release/pnetcdf-1.12.1.tar.gz

echo "##########build CMake##########"
tar xzf cmake-3.22.2.tar.gz
cd cmake-3.22.2
export CC=$c_compiler
export CXX=$cxx_compiler
export FC=$fortran_compiler
./configure --prefix=$installation_home/packages/cmake
make -j$procs && make install
cd ..

#echo "##########build mpich##########"
#tar xzf mpich-3.3.2.tar.gz
#cd mpich-3.3.2
#export CC=$c_compiler
#export CXX=$cxx_compiler
#export FC=$fortran_compiler
#./configure --prefix=$installation_home/packages/mpich
#make -j$procs && make install
#cd ..

echo "##########build parallel szip##########"
tar xzf szip-2.1.1.tar.gz
cd szip-2.1.1
export CC=$c_compiler
export CXX=$cxx_compiler
export FC=$fortran_compiler
./configure --prefix=$installation_home/packages/szip
make -j$procs && make install
cd ..

echo "##########build serial hdf5##########"
# build serial hdf5
tar xzf hdf5-1.10.6.tar.gz
cd hdf5-1.10.6
export CC=$c_compiler
export CXX=$cxx_compiler
export F77=$fortran_compiler
export FC=$fortran_compiler
./configure --prefix=$installation_home/packages/hdf5-serial --with-szlib=$installation_home/packages/szip --disable-parallel
make -j$procs && make install
cd ..

echo "##########build parallel hdf5##########"
rm -rf hdf5-1.10.6
tar xzf hdf5-1.10.6.tar.gz
cd hdf5-1.10.6
export CC=$mpic_compiler
export MPICC=$mpic_compiler
export CXX=$mpicxx_compiler
export MPICXX=$mpicxx_compiler
export FC=$mpifortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
./configure --prefix=$installation_home/packages/hdf5-parallel --with-szlib=$installation_home/packages/szip --enable-parallel
make -j$procs && make install
cd ..

echo "##########build serial netcdf-c##########"
# build serial netcdf
tar xzf v4.7.3.tar.gz
cd netcdf-c-4.7.3
export CC=$c_compiler
export MPICC=$mpic_compiler
export CXX=$cxx_compiler
export MPICXX=$mpicxx_compiler
export F77=$fortran_compiler
export MPIF77=$mpifortran_compiler
export FC=$fortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-I${installation_home}/packages/hdf5-serial/include -I${installation_home}/packages/szip/include"
export CXXFLAGS="-I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export LDFLAGS="-L$installation_home/packages/hdf5-serial/lib -L$installation_home/packages/szip/lib"
#export CFLAGS='-I/home/omarahme/original/packages/hdf5-serial/include -I/home/omarahme/original/packages/szip/include'
#export CXXFLAGS='-I/home/omarahme/original/packages/hdf5-serial/include -I/home/omarahme/original/packages/szip/include'
#export LDFLAGS='-L/home/omarahme/original/packages/hdf5-serial/lib -L/home/omarahme/original/packages/szip/lib'
./configure --prefix=$installation_home/packages/netcdf-serial --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-serial --with-szlib=$installation_home/packages/szip
make -j$procs && make install
cd ..

echo "##########build serial netcdf-cxx##########"
tar xzf v4.3.1.tar.gz
cd netcdf-cxx4-4.3.1
export CC=$c_compiler
export MPICC=$mpic_compiler
export CXX=$cxx_compiler
export MPICXX=$mpicxx_compiler
export F77=$fortran_compiler
export MPIF77=$mpifortran_compiler
export FC=$fortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export CXXFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export LDFLAGS="-L$installation_home/packages/netcdf-serial/lib -L$installation_home/packages/hdf5-serial/lib -L$installation_home/packages/szip/lib"
export LIBS='-lnetcdf -lhdf5_hl -lhdf5 -ldl -lsz -lm -lz'
./configure --prefix=$installation_home/packages/netcdf-serial --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-serial --with-szlib=$installation_home/packages/szip
make -j$procs && make install
cd ..

echo "##########build serial netcdf-fortran##########"
tar xzf v4.4.4.tar.gz
cd netcdf-fortran-4.4.4
export CC=$c_compiler
export MPICC=$mpic_compiler
export CXX=$cxx_compiler
export MPICXX=$mpicxx_compiler
export F77=$fortran_compiler
export MPIF77=$mpifortran_compiler
export FC=$fortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export CXXFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export FFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export FCFLAGS="-static -I$installation_home/packages/netcdf-serial/include -I$installation_home/packages/hdf5-serial/include -I$installation_home/packages/szip/include"
export LDFLAGS="-L$installation_home/packages/netcdf-serial/lib -L$installation_home/packages/hdf5-serial/lib -L$installation_home/packages/szip/lib"
export LIBS='-lnetcdf -lhdf5_hl -lhdf5 -ldl -lsz -lm -lz'
./configure --prefix=$installation_home/packages/netcdf-serial  --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-serial --with-szlib=$installation_home/packages/szip --enable-valgrind-tests --enable-serial-tests --enable-extra-tests --enable-extra-example-tests
make -j$procs && make install
cd ..

echo "##########build parallel netcdf-c##########"
# build parallel netcdf
rm -rf netcdf-c-4.7.3
tar xzf v4.7.3.tar.gz
cd netcdf-c-4.7.3
export CC=$mpic_compiler
export MPICC=$mpic_compiler
export CXX=$mpicxx_compiler
export MPICXX=$mpicxx_compiler
unset F77
unset MPIF77
export FC=$mpifortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export CXXFLAGS="-I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
unset FFLAGS
unset FCFLAGS
export LDFLAGS="-L$installation_home/packages/hdf5-parallel/lib -L$installation_home/packages/szip/lib"
unset LIBS
./configure --prefix=$installation_home/packages/netcdf-parallel --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-parallel --with-szlib=$installation_home/packages/szip
make -j$procs && make install
cd ..

echo "##########build parallel netcdf-cxx##########"
rm -rf netcdf-cxx4-4.3.1
tar xzf v4.3.1.tar.gz
cd netcdf-cxx4-4.3.1
export CC=$mpic_compiler
export MPICC=$mpic_compiler
export CXX=$mpicxx_compiler
export MPICXX=$mpicxx_compiler
export FC=$mpifortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export CXXFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export LDFLAGS="-L$installation_home/packages/netcdf-parallel/lib -L$installation_home/packages/hdf5-parallel/lib -L$installation_home/packages/szip/lib"
export LIBS='-lnetcdf -lhdf5_hl -lhdf5 -ldl -lsz -lm -lz'
./configure --prefix=$installation_home/packages/netcdf-parallel --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-parallel --with-szlib=$installation_home/packages/szip
make -j 4 && make install
cd ..

echo "##########build parallel netcdf-fortran##########"
rm -rf netcdf-fortran-4.4.4
tar xzf v4.4.4.tar.gz 
cd netcdf-fortran-4.4.4
export CC=$mpic_compiler
export MPICC=$mpic_compiler
export CXX=$mpicxx_compiler
export MPICXX=$mpicxx_compiler
export F77=$mpifortran_compiler
export MPIF77=$mpifortran_compiler
export FC=$mpifortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
export CFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export CXXFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export FFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export FCFLAGS="-static -I$installation_home/packages/netcdf-parallel/include -I$installation_home/packages/hdf5-parallel/include -I$installation_home/packages/szip/include"
export LDFLAGS="-L$installation_home/packages/netcdf-parallel/lib -L$installation_home/packages/hdf5-parallel/lib -L$installation_home/packages/szip/lib"
export LIBS='-lnetcdf -lhdf5_hl -lhdf5 -ldl -lsz -lm -lz'
./configure --prefix=$installation_home/packages/netcdf-parallel  --disable-dap --enable-netcdf-4 --disable-shared --enable-static --with-hdf5=$installation_home/packages/hdf5-parallel --with-szlib=$installation_home/packages/szip --enable-valgrind-tests --enable-parallel-tests --enable-extra-tests --enable-extra-example-tests
make -j$procs && make install
cd ..


echo "##########build pnetcdf##########"
# build pnetcdf
tar xzf pnetcdf-1.12.1.tar.gz
cd pnetcdf-1.12.1
export CC=$mpic_compiler
export MPICC=$mpic_compiler
export CXX=$mpicxx_compiler
export MPICXX=$mpicxx_compiler
unset F77
unset MPIF77
export FC=$mpifortran_compiler
export MPIF90=$mpifortran_compiler
export MPIFC=$mpifortran_compiler
unset CFLAGS
unset CXXFLAGS
unset FFLAGS
unset FCFLAGS
unset LDFLAGS
unset LIBS
./configure --prefix=$installation_home/packages/pnetcdf --disable-in-place-swap
make -j$procs && make install
