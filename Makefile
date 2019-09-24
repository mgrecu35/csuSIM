#-------------------------------------------------------------------------
#
#        Makefile for TB simulator from a combination of MIRS and CMB 
#
#-------------------------------------------------------------------------

FC      = ifort

FFLAGS  = -O2 -DLANGUAGE_FORTRAN -fPIC

INCFLAGS = -Ifastem/fastem_inc
LFLAGS = -Lfastem/fastem_lib -lm -lfastem

EXE     = GPM_TbsimulatorV7

#--------------------------------------------------------------------------
OBJS=\
	GPM_TbsimulatorV7_variables.o \
	GPM_TbsimulatorV7_IO.o \
	MONOrtm.o \
	Fastem.o \
	GPM_TbsimulatorV7_rad_procedures.o \
	GPM_TbsimulatorV7_mie_code.o \
	GPM_TbsimulatorV7_procedures.o \
	GPM_TbsimulatorV7.o readTables_nonsph.o retTablesInt.half.o\
	getScattProp.o bisection.o

OBJS2=\
	GPM_TbsimulatorV7_variables.o \
	GPM_TbsimulatorV7_IO.o \
	MONOrtm.o \
	GPM_TbsimulatorV7_rad_procedures.o \
	GPM_TbsimulatorV7_mie_code.o \
	GPM_TbsimulatorV7_procedures.o\
	getCSUvars.o readTables_nonsph.o retTablesInt.half.o\
	getScattProp.o bisection.o
			
#---------------------------------------------------------------------------
GPM_TbsimulatorV7: $(OBJS2)
	$(FC) $(FFLAGS) $(OBJS) $(LFLAGS) -o $(EXE) 
	ifort -openmp -shared -o \
		csuSIM.so $(OBJS2) -lm 
	#rm -f *.o *.mod

.SUFFIXES: .f90 .o	
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
readTables_nonsph.o: cmbCodes/readTables_nonsph.f90
	$(FC) -c $(FFLAGS)  cmbCodes/readTables_nonsph.f90

retTablesInt.half.o: cmbCodes/retTablesInt.half.f90
	$(FC) -c $(FFLAGS)  cmbCodes/retTablesInt.half.f90

getScattProp.o: cmbCodes/getScattProp.f90
	$(FC) -c $(FFLAGS)  cmbCodes/getScattProp.f90

bisection.o: cmbCodes/bisection.f90
	$(FC) -c $(FFLAGS)  cmbCodes/bisection.f90