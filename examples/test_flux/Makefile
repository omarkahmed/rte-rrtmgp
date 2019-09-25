#
# Here set variables RRTMGP_DIR, NCHOME, NFHOME, TIME_DIR (for GPTL)
#   or have those variables set in the environment
#
-include Makefile.libs
-include $(RRTMGP_DIR)/Makefile.conf
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_DIR)
TIME_DIR = /OSS_Scratch/RRTMGP/gptl

export RRTMGP_ROOT = /OSS_Scratch/RRTMGP/rte-rrtmgp/
export RRTMGP_BUILD = $(RRTMGP_ROOT)/build

VPATH = ../:./:$(RRTMGP_BUILD):$(RRTMGP_ROOT)/examples/:$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/extensions/cloud_optics

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib -L${HDF5}/lib
LIBS      += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5

#
# Setting macro TIMING=yes uses routines from the General Purpose Timing Library
#  https://jmrosinski.github.io/GPTL/
#
TIMING=yes
# Compiler specific - e.g. turn off for pgfortran, set to -cpp for gfortran
#FCFLAGS += -fpp
ifeq ($(TIMING),yes)
	#
	# Timing library
	#
	FCINCLUDE += -I$(TIME_DIR)/include
	# Compiler specific
	FCFLAGS += -DUSE_TIMING
	LDFLAGS   += -L$(TIME_DIR)/lib
	LIBS      += -lgptl
endif

DSRCPATH := ${PWD}/disort_lw
DLWSRC = \
	$(DSRCPATH)/rtregcld_ip.f \
	$(DSRCPATH)/rtregcld.f \
	$(DSRCPATH)/RDI1MACH.f \
	$(DSRCPATH)/LINPAK.f \
	$(DSRCPATH)/ErrPack.f \
	$(DSRCPATH)/rtrdis.f \
	$(DSRCPATH)/disort_lw.f
	
DLWSRCOBJ = $(patsubst %.f,%.o,$(DLWSRC))

# Compilation rules
%.o: %.f
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $< -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Ancillary codes
#
ADDITIONS = \
	${DLWSRCOBJ} \
	mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o mo_test_files_io.o mo_cloud_optics.o mo_load_cloud_coefficients.o  


all: test_lw

test_lw:  test_lw.o   $(ADDITIONS) $(RRTMGP_DIR)/librte.a $(RRTMGP_DIR)/librrtmgp.a

test_lw.o:  test_lw.F90   $(ADDITIONS) $(RRTMGP_DIR)/librte.a $(RRTMGP_DIR)/librrtmgp.a

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

mo_load_cloud_coefficients.o: mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.F90

disort_lw.o: ${DLWSRCOBJ}	

clean:
	-find . -iname '*.mod' -delete
	-find . -iname '*.o' -delete
	-rm -f test_lw 
