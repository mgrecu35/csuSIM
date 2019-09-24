#-------------------------------------------------------------------------
#
#            Makefile for stand alone version of FASTEM4/5/6
#
#-------------------------------------------------------------------------

CODE = fortcode
FC      = ifort

CFLAGS= -g -assume byterecl -O3

LIB =  -lm

AR = ar rc

FLAGS = $(LIB)

#--------------------------------------------------------------------------
OBJS=	Type_Kinds.o \
	Compare_Float_Numbers.o \
	Endian_Utility.o \
	File_Utility.o \
	Message_Handler.o \
	Binary_File_Utility.o \
	Date_Utility.o \
	Search_Utility.o \
	CRTM_Parameters.o \
	CRTM_Interpolation.o \
	Hyperbolic_Step.o \
	Slope_Variance.o \
	Fundamental_Constants.o \
	Fresnel.o \
	Liu.o \
	SensorInfo_Parameters.o \
	Sort_Utility.o \
	Subset_Define.o \
	ACCoeff_Define.o \
	NLTECoeff_Define.o \
	SpcCoeff_Define.o \
	ACCoeff_Binary_IO.o \
	NLTECoeff_Binary_IO.o \
	SpcCoeff_Binary_IO.o \
	Guillou.o \
	Ellison.o \
	CRTM_Fastem1.o \
	CRTM_SensorData_Define.o \
	CRTM_Surface_Define.o \
	CRTM_Geometry_Define.o \
	CRTM_GeometryInfo_Define.o \
	CRTM_GeometryInfo.o \
	CRTM_LowFrequency_MWSSEM.o \
	CRTM_SpcCoeff.o \
	CRTM_SfcOptics_Define.o \
	FitCoeff_Define.o \
	Foam_Utility_Module.o \
	Small_Scale_Correction_Module.o \
	Large_Scale_Correction_Module.o \
	Reflection_Correction_Module.o \
	MWwaterLUT_Define.o \
	MWwaterCoeff_Define.o \
	Azimuth_Emissivity_F6_Module.o \
	Azimuth_Emissivity_Module.o \
	CRTM_FastemX.o \
	CRTM_MWwaterCoeff.o \
	CRTM_MW_Water_SfcOptics.o

#---------------------------------------------------------------------------
libfastem.a: $(OBJS)
	$(AR) $@ $(OBJS)

.SUFFIXES: .f90 .f .o

.f90.o:
	$(FC) -c $(CFLAGS) $<
.f.o:
	$(FC) -c $(CFLAGS) $<

# clean up files

clean :
	-rm *.o *.mod
	-rm libfastem.a

