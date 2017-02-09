# Marco Santia
#
 
F90 = gfortran
FFLAGS = -c -o $@

poisson.x: math.o linalg.o globals.o poisson.o
	${F90} -o poisson.x math.o linalg.o globals.o poisson.o 

math.o: math.f90 
	${F90} ${FFLAGS} math.f90

linalg.o: linalg.f90 
	${F90} ${FFLAGS} linalg.f90

globals.o: globals.f90 
	${F90} ${FFLAGS} globals.f90

poisson.o: poisson.f90 
	${F90} ${FFLAGS} poisson.f90

clean:
	rm -f poisson.x math.o linalg.o globals.o poisson.o
	rm -f *.mod
	@echo "all cleaned up!"


