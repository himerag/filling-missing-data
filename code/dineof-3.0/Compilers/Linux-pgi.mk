#
# Include file for PGI Fortran compiler on Linux
#


F90C := pgf90
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -byteswapio
endif  

ifdef DEBUG
  F90FLAGS += -g -C
else
# -Mipa=fast has some troubles on some systems (mine)
#  F90FLAGS += -u -fastsse -Mipa=fast
# now -u -fastsse in both
  F90FLAGS += -O3
  LDFLAGS += 
endif

ifdef DOUBLEPREC
  F90FLAGS += -r8 -DDOUBLE_PRECISION
endif

ifdef STATIC
  F90FLAGS += -Bstatic
  LDFLAGS += -Bstatic
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

#NETCDF_LIB ?= $(HOME)/netcdf/lib
NETCDF_LIB ?= 
#==========================================================================
# This directory should contains a file usually called netcdf.inc and 
# netcdf.mod
#==========================================================================

#NETCDF_INC ?= $(HOME)/netcdf/include
NETCDF_INC ?= 
#==========================================================================
# All include options
#==========================================================================

INCLUDES = -I$(NETCDF_INC)

#==========================================================================
# All library link options
#==========================================================================

LIBRARIES = -L$(ARPACK_LIB) -larpack -llapack -lblas -L$(NETCDF_LIB) -lnetcdf
