SUBROUTINES = distances.f90 0thEnergy.f90
READ = read_parameters.f90 read_coordinates.f90 read_basis.f90

all:
	gfortran GFN1-xTB.f90 -o xTB.x ${SUBROUTINES} ${READ}

clean:
	rm -f xTB.x
