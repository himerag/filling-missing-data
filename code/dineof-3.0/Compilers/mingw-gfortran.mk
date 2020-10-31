#==========================================================================
#
# Include file for gfortran Fortran compiler on Windows (mingw)
# For more information about mingw:
# http://www.mingw.org
# http://gcc.gnu.org/wiki/GFortranBinariesWindows
#
#==========================================================================


F90C := gfortran

F90FLAGS :=

LD := $(F90C)
LDFLAGS := 

# necessary for gfortran 4.1, but default since 4.2
F90FLAGS += -frecord-marker=4


ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian
endif  

ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif

ifdef STATIC
  LDFLAGS += -static
endif

ifdef DOUBLEPREC
  F90FLAGS += -fdefault-real-8 -DDOUBLE_PRECISION
endif

#==========================================================================
# Library locations
#==========================================================================

#==========================================================================
# Locations of the ARPACK library
# This directory should contains a file usually called libarpack.a
#==========================================================================

ARPACK_LIB=

#==========================================================================
# Locations of the NetCDF library
# This directory should contains a file usually called libnetcdf.a
#==========================================================================

NETCDF_LIB=

#==========================================================================
# This directory should contains a file usually called netcdf.inc and 
# netcdf.mod
#==========================================================================

NETCDF_INC=

#==========================================================================
# All include options
#==========================================================================

INCLUDES = -I$(NETCDF_INC)

#==========================================================================
# All library link options
#==========================================================================

LIBRARIES = -L$(ARPACK_LIB) -larpack -L$(NETCDF_LIB) -lnetcdf
