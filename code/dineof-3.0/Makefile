

#==========================================================================
#
# Implemented combinations of OS and FORT are:
#
#  +--------+-----------+
#  | OS     | FORT      |
#  |--------+-----------+
#  | Linux  | ifort     |
#  | Linux  | pgi       |
#  | Linux  | g95       |
#  | Linux  | gfortran  |
#  | IRIX   | f90       |
#  +--------+-----------+
#
# For ports to other achictectures, create a corresponding file
# in the "Compilers" directory. Contributions are most welcome!.
#
#==========================================================================

OS ?= Linux
FORT ?= gfortran

#==========================================================================
# Start of user-defined options. Modify macro variables: 'on' is TRUE while
# blank is FALSE.
#==========================================================================
#
#  Activate debugging compiler options:

DEBUG ?= 

# DINEOF version
VERSION ?= 3.0

#  Activate static linking

STATIC ?= 

#==========================================================================
#  Machine format, can be either big_endian, little_endian or native machine
#  format. YOUR CHOICE MUST BE CONSITENT WITH THE WAY YOU CREATE THE INPUT
#  FILE (see also matlab routines gread.m and gwrite.m)
#
#  Note: with some compilers (e.g. IRIX-f90) some options are
#  not available
#==========================================================================

#FORMAT := little_endian
FORMAT := big_endian
#FORMAT := native

#==========================================================================
# PRECISION attributed for the internal data processing
#==========================================================================

DOUBLEPREC ?= ON

#==========================================================================
# Include compiler specific options
#==========================================================================

include Compilers/$(OS)-$(strip $(FORT)).mk

F90FLAGS += $(INCLUDES)

#==========================================================================
# All Source files
#==========================================================================

SOURCE = ReadMatrix.F90 ufileformat.F90 initfile.F90 stat.F90 norm.F90 dineof_utils.F90 \
  smeanToZero.F smeanByRow.F svariExp.F ssvd_lancz.F valsvd.F90 dineof.F90 writeSVD.F90 \
  writeMatrix.F90

#==========================================================================
# Object files for "dineof"
#==========================================================================


OBJECTS =  ReadMatrix.o ufileformat.o initfile.o stat.o norm.o dineof_utils.o \
  smeanToZero.o smeanByRow.o svariExp.o ssvd_lancz.o valsvd.o dineof.o writeSVD.o \
  writeMatrix.o 

#==========================================================================
# Object files for "crossval"
#==========================================================================


OBJECTSC =  ReadMatrix.o ufileformat.o initfile.o stat.o norm.o crossval.o

#==========================================================================
# Declare .F90 a valid suffix
#==========================================================================

.SUFFIXES: .F90 .F

#==========================================================================
# How to compile F77 and F90 programs?
#==========================================================================

%.o: %.F
	$(F90C) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

all: crossval dineof

dineof: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBRARIES)

crossval: $(OBJECTSC)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTSC) $(LIBRARIES)

clean:
	rm -f $(OBJECTS) $(OBJECTSC) *.mod 

distclean:
	rm -f  $(OBJECTS) $(OBJECTSC) *.mod dineof crossval


release:
	rm -rf /tmp/dineof-$(VERSION);\
	svn export https://dineof.googlecode.com/svn/trunk /tmp/dineof-$(VERSION);\
	tar -C /tmp/ -zcvf dineof-$(VERSION).tar.gz dineof-$(VERSION);\
	rm -rf /tmp/dineof-$(VERSION)

tarfile:
	if [ -e dineof-$(VERSION).tar.gz ]; then  rm -i dineof-$(VERSION).tar.gz; fi
	tar -C ../ --exclude-vcs -zcvf dineof-$(VERSION).tar.gz dineof

print:
	echo $(OBJECTS) 

#==========================================================================
# Dependencies
#==========================================================================

ufileformat.o: ppdef.h
initfile.o: ppdef.h
crossval.o: ufileformat.o initfile.o
ReadMatrix.o: ufileformat.o
writeMatrix.o: ufileformat.o
writeSVD.o:  ufileformat.o
dineof.o: ReadMatrix.h ufileformat.o initfile.o 
ssvd_lancz.f: includes/debug.h 
stat.o: ufileformat.o
siterativeEof.o:  ufileformat.o
ssvd_lancz.o: ufileformat.o dineof_utils.o
dineof_utils.o: ufileformat.o
