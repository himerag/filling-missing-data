#==========================================================================
#
# Include file for PGI Fortran compiler on Linux
#
#==========================================================================


F90C := ifort
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -convert big_endian
endif  

ifdef DEBUG
  F90FLAGS += -g
else
#  F90FLAGS += -fast 
endif

ifdef DOUBLEPREC
  F90FLAGS += -r8 -DDOUBLE_PRECISION
endif


#==========================================================================
# Library locations
#==========================================================================

#==========================================================================
# Locations of the ARPACK library
# This directory should contains a file usually called libarpack.a
#==========================================================================

ARPACK_LIB ?=

#==========================================================================
# Locations of the NetCDF library
# This directory should contains a file usually called libnetcdf.a
#==========================================================================

NETCDF_LIB ?=

#==========================================================================
# This directory should contains a file usually called netcdf.inc and 
# netcdf.mod
#==========================================================================

NETCDF_INC ?=

#==========================================================================
# All include options
#==========================================================================

INCLUDES = -I$(NETCDF_INC)

#==========================================================================
# All library link options
#==========================================================================

LIBRARIES = -L$(ARPACK_LIB) -larpack_ifort -llapack -lblas -L$(NETCDF_LIB) -lnetcdf
