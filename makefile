#code to deal with the yosemette compile issue

CC=gfortran
CFLAGS=-cpp
BASE=$(PWD)
BIN=$(BASE)/bin
STAR=$(BIN)/stellarstructure
ADIABAT=$(BIN)/stellar_adiabat
SHOOT=$(BIN)/shooting_solution
SRC=$(BASE)/src
SOLVERSTEST=$(BIN)/solvers_unit_test

all: $(STAR) $(ADIABAT) $(SHOOT)

$(STAR): $(SRC)/mod_solvers.o $(SRC)/jacobian_solution.o $(SRC)/read_io.o
	mkdir -p bin
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) jacobian_solution.o mod_solvers.o read_io.o mod_shared.o -o $(STAR)

$(ADIABAT): $(SRC)/mod_solvers.o $(SRC)/adiabatic_solution.o $(SRC)/read_io.o
	mkdir -p bin
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) adiabatic_solution.o mod_solvers.o read_io.o mod_shared.o -o $(ADIABAT)

$(SHOOT): $(SRC)/mod_solvers.o $(SRC)/shooting_solution.o $(SRC)/read_io.o
	mkdir -p bin
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) shooting_solution.o mod_solvers.o read_io.o mod_shared.o -o $(SHOOT)

$(SRC)/mod_solvers.o: $(SRC)/mod_solvers.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) $(CFLAGS) -c mod_solvers.f90

$(SRC)/shooting_solution.o: $(SRC)/shooting_solution.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) -c shooting_solution.f90 mod_shared.o

$(SRC)/jacobian_solution.o: $(SRC)/jacobian_solution.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) -c jacobian_solution.f90 mod_shared.o

$(SRC)/adiabatic_solution.o: $(SRC)/adiabatic_solution.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) -c adiabatic_solution.f90 mod_shared.o

$(SRC)/read_io.o: $(SRC)/read_io.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) $(CFLAGS) -I$(SRC) -c read_io.f90 mod_shared.o

$(SRC)/mod_shared.o: $(SRC)/mod_shared.f90 
	cd $(SRC); $(CC) $(CFLAGS) -c mod_shared.f90

#split this up into seperate unittests when necisary
test: $(SRC)/mod_solvers.f90
	cd $(SRC); $(CC) $(CFLAGS) -fbounds-check -Dunittest mod_solvers.f90 -o $(SOLVERSTEST)
	$(SOLVERSTEST)

clean: 
	rm src/*.o
	rm src/*.mod
