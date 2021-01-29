#include <netcdf>
#include "YAKL.h"

using namespace yakl;
namespace rrtmgp {

    class MySimpleNetCDF {


        protected:

            int ncid;

        public:

            // Constructor
            MySimpleNetCDF() {};

            // Destructor
            ~MySimpleNetCDF() {
                //close();
            };

            void close() {
                handle_error(nc_close(ncid));
            }

            void create(std::string filename, int mode=NC_CLOBBER) {
                handle_error(nc_create(filename.c_str(), mode, &ncid));
            };

            void open(std::string filename, int mode=NC_NOWRITE) {
                handle_error(nc_open(filename.c_str(), mode, &ncid));
            };

            void open(char *filename) {
                handle_error(nc_open(filename, NC_NOWRITE, &ncid));
            }

            // NetCDF routines return an integer error code. Define a function
            // here to abort program execution and throw an error code if we
            // encounter a non-zero NetCDF return code. We will wrap our
            // NetCDF calls with this function to handle these errors in a
            // consistent way
            void handle_error(int err) {
                if (err) {
                    std::cout << "ERROR: " << nc_strerror(err) << std::endl;
                    abort();
                }
            }

            void handle_error(int err, const char *file, int line) {
                if (err) {
                    std::cout << "ERROR: " << nc_strerror(err) << " at line " << line << " in " << file << std::endl;
                    abort();
                }
            }

            // Read a netCDF array to a YAKL array
            template <class T, int rank, int myMem, int myStyle> void read(Array<T,rank,myMem,myStyle> &arr, std::string varName) {

                // Get variable ID
                int varid;
                handle_error(nc_inq_varid(ncid, varName.c_str(), &varid));

                // Get variable dimension sizes
                int ndims;
                int dimids[NC_MAX_VAR_DIMS];
                nc_type vtype;
                handle_error(nc_inq_var(ncid, varid, NULL, &vtype, &ndims, dimids, NULL));
                std::vector<int> dimSizes(ndims);
                size_t dimsize;
                for (int i = 0; i < ndims; i++) {
                    handle_error(nc_inq_dimlen(ncid, dimids[i], &dimsize)); 
                    dimSizes[i] = dimsize;
                }

                // If style is fortran, we need to reverse array dims
                if (myStyle == styleFortran) {
                    std::reverse(dimSizes.begin(), dimSizes.end());
                }

                // Allocate (or reshape) the yakl array
                arr = Array<T,rank,myMem,myStyle>(varName.c_str(),dimSizes);

                // Read variable data
                if (myMem == memDevice) {
                    auto arrHost = arr.createHostCopy();
                    if (std::is_same<T,bool>::value) {
                        // Create boolean array from integer arrays
                        Array<int,rank,memHost,myStyle> tmp("tmp",dimSizes);
                        handle_error(nc_get_var(ncid, varid, tmp.data()));
                        for (int i=0; i < arr.totElems(); i++) { arrHost.myData[i] = tmp.myData[i] == 1; }
                    } else {
                        // Need to be careful with floats; nc_get_var is overloaded on type, but we need
                        // to make sure we read floats from file with the float procedure, and doubles
                        // with that for doubles. The danger is if the user passes a yakl array here
                        // with type double, but tries to read type float from file.
                        // TODO: why does the YAKL implementation for this work fine, but this version
                        // calling nc_get_var directly does not?
                        if (vtype == NC_FLOAT) {
                            Array<float,rank,memHost,myStyle> tmp("tmp",dimSizes);
                            handle_error(nc_get_var(ncid, varid, tmp.data()));
                            for (int i=0; i < arr.totElems(); i++) { arrHost.myData[i] = tmp.myData[i]; }
                        } else {
                            handle_error(nc_get_var(ncid, varid, arrHost.data()));
                        }
                    }
                    arrHost.deep_copy_to(arr);
                } else {
                    if (std::is_same<T,bool>::value) {
                        // Create boolean array from integer arrays
                        Array<int,rank,memHost,myStyle> tmp("tmp",dimSizes);
                        handle_error(nc_get_var(ncid, varid, tmp.data()));
                        for (int i=0; i < arr.totElems(); i++) { arr.myData[i] = tmp.myData[i] == 1; }
                    } else {
                        if (vtype == NC_FLOAT) {
                            Array<float,rank,memHost,myStyle> tmp("tmp",dimSizes);
                            handle_error(nc_get_var(ncid, varid, tmp.data()));
                            for (int i=0; i < arr.totElems(); i++) { arr.myData[i] = tmp.myData[i]; }
                        } else {
                            handle_error(nc_get_var(ncid, varid, arr.data()));
                        }
                    }
                }

            }

            // Read a scalar type
            template <class T> void read(T &arr , std::string varName) {
                // Get variable ID
                int varid;
                handle_error(nc_inq_varid(ncid, varName.c_str(), &varid));

                // Read data
                handle_error(nc_get_var(ncid, varid, &arr));
            }

            // Check if variable exists in file
            bool varExists (std::string varName) {
                int varid;
                int ncerr = nc_inq_varid(ncid, varName.c_str(), &varid);
                if (ncerr == 0) {
                    return true;
                } else {
                    return false;
                }
            }

            bool dimExists (std::string dimName) {
                int dimid;
                int ncerr = nc_inq_dimid(ncid, dimName.c_str(), &dimid);
                if (ncerr == 0) {
                    return true;
                } else {
                    return false;
                }
            }

            void addDim(std::string dimName, int dimSize) {
                // Put file into define mode
                int ncerr = nc_redef(ncid);
                if ((ncerr != NC_NOERR) and (ncerr != NC_EINDEFINE)) {
                    handle_error(ncerr);
                }

                // Define dimension
                int dimid;
                handle_error(nc_def_dim(ncid, dimName.c_str(), dimSize, &dimid));

                // End define mode
                handle_error(nc_enddef(ncid));
            }

            template <class T, int rank, int myMem, int myStyle> 
            void write(Array<T,rank,myMem,myStyle> const &arr, std::string varName, std::vector<std::string> dimNames) {

                // Make sure length of dimension names is equal to rank of array
                if (rank != dimNames.size()) { yakl_throw("dimNames.size() != Array rank"); }

                // Define dimensions if they do not exist
                for (int i = 0; i < dimNames.size(); i++) {
                    if (dimExists(dimNames[i])) {
                        // check that size is correct
                    } else {
                        addDim(dimNames[i], arr.dimension[i]);
                    }
                }

                // Write data to file
                addVar(arr, varName);
            }
    };

}
