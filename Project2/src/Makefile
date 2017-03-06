# Marco Santia
 
F90 = gfortran
FFLAGS = -O3 -pg -c -o $@

solver.x: linalg.o input.o solver.o
	${F90} -o solver.x linalg.o input.o solver.o 

linalg.o: linalg.f90 
	${F90} ${FFLAGS} linalg.f90

input.o: input.f90
	${F90} ${FFLAGS} input.f90

solver.o: solver.f90 
	${F90} ${FFLAGS} solver.f90

clean:
	rm -f solver.x linalg.o input.o solver.o
	rm -f *.mod
	@echo "all cleaned up!"


