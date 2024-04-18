maindir = ../
bin = bin/

SUBROUTINES = distances.f90 0thEnergy.f90 print_matrix.f90 inv_sqrt_S.f90 CK_matrix.f90 SCF.f90
READ = read_parameters.f90 read_coordinates.f90 read_basis.f90
LAPACK = lapack_routines.f90

GF = gfortran
LFLAGS = -llapack -lblas
DFLAGS = -g -Wall

program:
	$(GF) GFN1-xTB.f90 -o $(maindir)$(bin)xTB.x $(SUBROUTINES) $(READ) $(LAPACK) $(LFLAGS)

debug:
	$(GF) $(DFLAGS) GFN1-xTB.f90 -o $(maindir)$(bin)debug.x $(SUBROUTINES) $(READ) $(LAPACK) $(LFLAGS)

clean:
	rm -f $(bin)xTB.x
	rm -f $(bin)debug.x
