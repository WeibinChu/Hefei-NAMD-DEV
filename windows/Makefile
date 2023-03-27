#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
# FC       = gfortran
FC       = x86_64-w64-mingw32-gfortran

FFLAGS   = -g -O3 -fbacktrace -fbounds-check -fopenmp # -Wall -Wextra -Wconversion
MAKE     = make

BLAS     = windows/libblas.dll.a
LAPACK   = windows/liblapack.dll.a
MPI      = windows/libmsmpi.a
LLIBS    = $(BLAS) $(LAPACK) $(MPI)
DMPI     = -D "INT_PTR_KIND()=8" -fno-range-check -DENABLEMPI 
DMKL     = -DENABLEMKL

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC = prec.f90 utils.f90 parallel.f90 fileio.f90 couplings.f90 hamil.f90 \
   	TimeProp.f90 dish.f90 fssh.f90 main.f90


OBJ = $(SRC:.f90=.o)
EXE = hfnamd

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $< 

.SUFFIXES: .o .F90
.F90.o:
	$(FC) $(FFLAGS) -c $< $(DMPI) $(DMKL)


#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
tdm:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(LLIBS)

prepare:
	cp windows/mpi.mod .
	cp windows/mpi.F90 .
clean:
	rm -f *.mod
	rm -f *.o
	rm -f mpi.F90
tar:
	tar -czvf hfnamd.tgz *.f90 Makefile
tag:
	ctags *.f90
