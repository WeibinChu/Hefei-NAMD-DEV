#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
# FC= ifort -assume byterecl 
FC= gfortran
FC= x86_64-w64-mingw32-gfortran
# FFLAGS= -g -O3
FFLAGS= -g -O3 -fbacktrace -fbounds-check # -Wall -Wextra -Wconversion
# FFLAGS= -g -O3 -heap-arrays -check bounds -g -traceback
MAKE = make
LAPACK = libblas.dll.a liblapack.dll.a

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC = prec.f90 utils.f90 fileio.f90 couplings.f90 hamil.f90 \
   	TimeProp.f90 dish.f90 fssh.f90 main.f90


OBJ = $(SRC:.f90=.o)
EXE = hfnamd

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
tdm:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(SPGLIB) $(LAPACK) 

clean:
	rm -f *.mod
	rm -f $(OBJ)
tar:
	tar -czvf hfnamd.tgz *.f90 Makefile
tag:
	ctags *.f90
