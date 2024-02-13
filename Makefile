all:
	gfortran -c util_module.f90 chem_module.f90
	gfortran GFN1-xTB.f90 -o xTB.x util_module.o chem_module.o
	rm *.o
