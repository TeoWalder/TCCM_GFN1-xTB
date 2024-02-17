all:
	gfortran -c module_utils.f90 module_chem.f90
	gfortran GFN1-xTB.f90 -o xTB.x module_utils.o module_chem.o
	rm *.o
