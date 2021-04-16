#**********************************************************
#*                                                        *
#*      ANT.G-2.5.0 - Makefile                            *
#*                                                        *
#**********************************************************
#*                                                        *
#*  Copyright (c) by                                      *
#*                                                        *
#*  David Jacob (1)                                       *
#*  Juan Jose Palacios (2)                                *
#*                                                        *
#* (1) Theory Department                                  *
#*     Max-Planck-Institute for Microstructure Physics    *
#*     Halle, 06120 (GERMANY)                             *
#* (2) Departamento de Fisica de la Materia Condensada    *
#*     Universidad Autonoma de Madrid                     *      
#*     28049 Madrid (SPAIN)                               *
#**********************************************************

# set complier/linker and options
MAKE=make
F90=ifort
CC=icc
AR=ar
BIN=../bin

# MKL related libs
##MKLROOT= /opt/intel/composerxe-2011.2.137/mkl
MKLROOT=/software/intel/parallel_studio_xe_2016/compilers_and_libraries_2016.1.150/linux/mkl
MKLPATH= ${MKLROOT}/lib/intel64
MKLINCLUDE=${MKLROOT}/include
INTELLIBS= -Wl,--start-group ${MKLPATH}/libmkl_intel_ilp64.a ${MKLPATH}/libmkl_core.a ${MKLPATH}/libmkl_intel_thread.a -Wl,--end-group

include make.in

# use these sources: 
SRCS= preproc.F util.f90 parameters.f90 antcommon.f90 g09Common.f90 constants.f90 numeric.f90 cluster.f90 ortho.f90 correlation.f90 OneDLead.f90 BetheLattice.f90 SpinOrbit.f90 MolMod.f90 device.f90 ANT.f90

#SRCS= preproc.F util.f90 parameters.f90 antcommon.f90 g09Common.f90 constants.f90 lapack_blas.f numeric.f90 cluster.f90 ortho.f90 correlation.f90 OneDLead.f90 BetheLattice.f90 device.f90 ANT.f90

# use name of parent directory as library name
PARDIRNAME=$(shell dirname `pwd`)
LIBNAME=$(shell basename $(PARDIRNAME))

# objects
OBJS= $(SRCS:%.c=../obj/%.o) $(SRCS:%.f=../obj/%.o) $(SRCS:%.F=../obj/%.o) $(SRCS:%.f90=../obj/%.o) 

# ANT.G library to be generated
ANTLIB=../lib/$(LIBNAME).a
L502=$(BIN)/l502.exe 

# prepocessing 
%.f90: %.F90
	gau-cpp -DG09ROOT -DINTEL $< $@

# rule for object generation
../obj/%.o: %.c 
	$(CC) -c -o $@ $<
../obj/%.o: %.F
	$(F90) $(CFLAGS) -I$(MKLINCLUDE) $(PREPROC) -c -o $@ $<
../obj/%.o: %.f 
	$(F90) $(CFLAGS) -I$(MKLINCLUDE) -c -o $@ $<
../obj/%.o: %.f90 
	$(F90) $(CFLAGS) -I$(MKLINCLUDE) -c -o $@ $<

# targets: libraries and l502.exe
target: $(ANTLIB) $(L502) 

# rule for the library
$(ANTLIB): $(OBJS)
	echo $(OBJS)
	- rm  $(ANTLIB)
	$(AR) crv $(ANTLIB) $(OBJS)

# rule for l502.exe
$(L502): $(OBJS)
	$(F90) $(CFLAGS) -o $(L502) $(BIN)/ml502.o $(BIN)/l502.a $(BIN)/bdam1.o $(BIN)/caldsu.o -L$(MKLPATH) $(GAUSSLIBS) $(ANTLIB) $(INTELLIBS) $(LFLAGS) 

	chmod o-rwx $(L502)

# Cleaning up to compile from scratch afterwards
clean:
	@echo "Cleaning up..."
	-rm ./*.mod
	-rm ./*.f90
	-rm $(BIN)/l502.exe
	-rm ../obj/*.o
	-rm ../lib/*.a
