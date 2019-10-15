#*****************************************************
#  IMPACT-T Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran compiler and linker.

FC = mpifort
LINK = mpifort
FFLAGS = -O3

#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJECTS = \
	NumConst.o PhysConst.o Pgrid.o Data.o \
	Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o \
	BPM.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o Multipole.o ConstFoc.o \
	SolRF.o Sol.o Dipole.o EMfld.o EMfldCart.o EMfldCyl.o EMfldAna.o \
	RFQ.o Fieldmap.o \
	BeamLineElem.o CompDom.o BeamBunch.o Depositor.o Ranger.o Field.o Distribution.o \
	Input.o Output.o AccSimulator.o main.o	

#**************************************************************************
# Search path directives

vpath %.f90 Appl Contrl DataStruct Func

#**************************************************************************
# Define name and location of output executable

EXENAME = ImpactTexe
INSTALLDIR = ~/Simulations/Current

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

#**************************************************************************
# Rules for building EXENAME from OBJECTS

$(EXENAME):  $(OBJECTS)
	$(LINK) -o $(EXENAME) $(OBJECTS)

#************************************************************************
# Copy module code to subfolders for debugging

	cp  *.mod Contrl
	cp  *.mod Appl
	cp  *.mod Func
	cp  *.mod DataStruct

#**************************************************************************
# Dependencies of specific OBJECTS

Data.o : NumConst.o
NumConst.o :
Pgrid.o :
PhysConst.o :
Bessel.o :
FFT.o : Timer.o Transpose.o
Fldmger.o : Timer.o Pgrid.o
Ptclmger.o : Timer.o Pgrid.o
Timer.o :
Transpose.o : Timer.o
BeamBunch.o : Timer.o PhysConst.o BeamLineElem.o CompDom.o Pgrid.o \
							Fldmger.o Data.o
BeamLineElem.o : Quadrupole.o Multipole.o DriftTube.o CCL.o CCDTL.o DTL.o \
								 SC.o BPM.o ConstFoc.o SolRF.o Dipole.o EMfld.o EMfldAna.o \
								 EMfldCyl.o EMfldCart.o RFQ.o Fieldmap.o
BPM.o : PhysConst.o Data.o
CompDom.o : Timer.o Pgrid.o
ConstFoc.o : Timer.o Pgrid.o
Depositor.o : Pgrid.o CompDom.o BeamBunch.o PhysConst.o
Dipole.o : Timer.o Pgrid.o
Distribution.o : Pgrid.o Timer.o NumConst.o PhysConst.o BeamBunch.o CompDom.o
DriftTube.o : PhysConst.o Data.o
DTL.o : PhysConst.o Data.o
EMfld.o : PhysConst.o Data.o
EMfldAna.o : PhysConst.o Data.o
EMfldCyl.o : PhysConst.o Data.o
EMfldCart.o : PhysConst.o Data.o
Field.o : Timer.o CompDom.o Pgrid.o FFT.o PhysConst.o Transpose.o Fldmger.o
Fieldmap.o : PhysConst.o Data.o
Multipole.o : PhysConst.o Data.o
Quadrupole.o : PhysConst.o Data.o
Ranger.o : Pgrid.o
RFQ.o : PhysConst.o Data.o
SC.o : PhysConst.o Data.o
Sol.o : PhysConst.o Data.o
SolRF.o : PhysConst.o Data.o
AccSimulator.o : Pgrid.o CompDom.o Field.o BeamLineElem.o Ptclmger.o \
								 BeamBunch.o Timer.o Input.o Output.o Data.o PhysConst.o \
								 NumConst.o Distribution.o Ranger.o Depositor.o
Input.o :
Output.o : Timer.o BeamBunch.o PhysConst.o
main.o : AccSimulator.o

#***********************************************************************
install: $(EXENAME)
	cp $(EXENAME) $(INSTALLDIR)/
clean:
	-rm -f $(EXENAME) *.o *.mod */*.o */*.mod
