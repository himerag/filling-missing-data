#==========================================================================
#
# Include file for g95 Fortran compiler on Linux
#
#==========================================================================


F90C := g95
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fendian=big
else
  ifeq ($(FORMAT),little_endian)
    F90FLAGS += -fendian=little
  endif
endif  

ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif

ifdef STATIC
  LDFLAGS += -static
endif

#==========================================================================
# Library locations
#==========================================================================

#==========================================================================
# Locations of the ARPACK library
# This directory should contain a file usually called libarpack.a
#==========================================================================

ARPACK_LIB ?= $(HOME)/ARPACK

#==========================================================================
# Locations of the NetCDF library
# This directory should contains a file usually called libnetcdf.a
#==========================================================================

NETCDF_LIB ?= /usr/local/lib

#==========================================================================
# This directory should contains a file usually called netcdf.inc and 
# netcdf.mod
#==========================================================================

NETCDF_INC ?= /usr/local/include/

#==========================================================================
# All include options
#==========================================================================

INCLUDES = -I$(NETCDF_INC)

#==========================================================================
# All library link options
#==========================================================================

LIBRARIES = -L$(ARPACK_LIB) -larpack -llapack -lblas -L$(NETCDF_LIB) -lnetcdf
