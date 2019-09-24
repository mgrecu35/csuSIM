#-------------------------------------------------------------------------
#
#            Makefile for stand alone version of FASTEM4/5/6
#
#-------------------------------------------------------------------------

CODE = fortcode
FC      = ifort

CFLAGS= -g

#INCFLAGS = -I/home/berg/src/fastem/fastem_lib
INCFLAGS = -I/home/drandel/gpm/dbase/mkdbaseV6/processV6/fastem/fastem_inc
LFLAGS = -L/home/drandel/gpm/dbase/mkdbaseV6/processV6/fastem/fastem_lib -lm -lfastem

EXE_NAME     = cmpr_emis

#--------------------------------------------------------------------------
OBJS=	RSS_RTM_module.o \
	Fastem.o \
	cmpr_emis.o

#---------------------------------------------------------------------------
cmpr_emis:\
	$(OBJS)
	$(FC) $(OBJS) $(LFLAGS) -o $(EXE_NAME) -check all

	rm -f *.o *.mod

.SUFFIXES: .f90 .f .o

.f90.o:
	$(FC) -c $(CFLAGS) $(INCFLAGS) $<
.f.o:
	$(FC) -c $(CFLAGS) $(INCFLAGS) $<

