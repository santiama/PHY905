# Marco Santia
#
 
F90 = gfortran
FFLAGS = -c -o $@

solver.x: input.o math.o linalg.o HarmonicOscillator.o poisson.o solver.o
	${F90} -o solver.x input.o math.o linalg.o HarmonicOscillator.o poisson.o solver.o 

input.o: input.f90
	${F90} ${FFLAGS} input.f90

math.o: math.f90 
	${F90} ${FFLAGS} math.f90

linalg.o: linalg.f90 
	${F90} ${FFLAGS} linalg.f90

poisson.o: poisson.f90 
	${F90} ${FFLAGS} poisson.f90

HarmonicOscillator.o: HarmonicOscillator.f90 
	${F90} ${FFLAGS} HarmonicOscillator.f90

solver.o: solver.f90 
	${F90} ${FFLAGS} solver.f90

clean:
	rm -f solver.x input.o math.o linalg.o HarmonicOscillator.o poisson.o solver.o
	rm -f *.mod
	@echo "all cleaned up!"


