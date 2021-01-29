#include "simple_netcdf.hpp"
#include "const.h"
#include "YAKL_netcdf.h"

int main(int argc, char *argv[]) {
    yakl::init();
    if (argc < 2) {
        std::cout << "Need to specify filename." << std::endl;
        return 1;
    }
    auto filename = argv[1];

    // Create some dummy data
    int nx = 2;
    int ny = 3;
    real2d myarr("myarr", nx, ny);
    for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
            myarr(i, j) = i + (j - 1) * nx;
        }
    }
    float myscl = 5;

    // Create a dummy file and write data to dummy file
    std::vector<std::string> dimNames = {"nx", "ny"};
    simple_netcdf::SimpleNetCDF io;
    io.create("foo.nc", NC_CLOBBER);
    io.write(myarr, "myarr", dimNames);
    io.write(myscl, "myscl");
    io.close();

    // Try to read back data from dummy file and confirm it matches what we wrote
    real2d myarr2("myarr");
    float myscl2;
    io.open("foo.nc", NC_NOWRITE);
    io.read(myarr2, "myarr");
    io.read(myscl2, "myscl");
    io.close();

    for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
            if (myarr(i, j) != myarr2(i, j)) {
                std::cout << "Arrays not equivalent at " << i << ", " << j << "; " << myarr(i,j) << " != " << myarr2(i,j) << std::endl;
                return -1;
            }
        }
    }
    if (myscl2 != myscl) {
        std::cout << "Scalars are not equivalent: " << myscl2 << " != " << myscl << std::endl;
        return -1;
    }

    // Read something from RRTMGP file
    real1d solar_source;
    io.open(filename, NC_NOWRITE);
    io.read(solar_source, "solar_source");
    io.close();

    return 0;

}
