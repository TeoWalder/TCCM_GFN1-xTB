SUBROUTINES = distances.f90 0thEnergy.f90 print_matrix.f90 inv_sqrt_S.f90 CK_matrix.f90 SCF.f90
READ = read_parameters.f90 read_coordinates.f90 read_basis.f90
LAPACK = lapack_routines.f90 -llapack -lblas

GF = gfortran

all:
	${GF} GFN1-xTB.f90 -o xTB.x ${SUBROUTINES} ${READ} ${LAPACK}

clean:
	rm -f xTB.x
