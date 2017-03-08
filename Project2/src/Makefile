# Marco Santia
 
F90 = gfortran
FFLAGS = -pg -c -o $@

solver.x: linalg.o unit_tests.o input.o solver.o
	${F90} -pg -o solver.x linalg.o unit_tests.o input.o solver.o 

linalg.o: linalg.f90 
	${F90} ${FFLAGS} linalg.f90

unit_tests.o: unit_tests.f90
	${F90} ${FFLAGS} unit_tests.f90

input.o: input.f90
	${F90} ${FFLAGS} input.f90

solver.o: solver.f90 
	${F90} ${FFLAGS} solver.f90

clean:
	rm -f solver.x linalg.o unit_tests.o input.o solver.o
	rm -f *.mod
	@echo "all cleaned up!"


