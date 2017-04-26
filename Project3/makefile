F90 = gfortran
FFLAGS = 
OBJ = bodies_class.o program3.o globals.o
MODULES = bodies_class.mod globals.mod
BIN = project3 

.PHONY: all clean

all: $(BIN)

$(BIN): $(OBJ)
	$(F90) -o $(BIN) $(OBJ)

globals.o globals.mod: globals.f90
	$(F90) -c $(FFLAGS) globals.f90

bodies_class.o bodies_class.mod: bodies_class.f90 globals.mod globals.o 
	$(F90) -c $(FFLAGS) bodies_class.f90

program3.o: program3.f90 globals.mod globals.o bodies_class.mod bodies_class.o
	$(F90) -c $(FFLAGS) program3.f90

clean:
	rm -f $(OBJ) $(MODULES)
