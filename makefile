CC=gfortran
BASE=$(PWD)
STAR=$(BASE)/bin/stellarstructure
SRC=$(BASE)/src

all: $(STAR)

$(STAR): $(SRC)/solvers.o $(SRC)/shooting_solution.o $(SRC)/read_io.o
	mkdir -p bin
	cd $(SRC); $(CC) -I$(SRC) shooting_solution.o mod_solvers.o read_io.o mod_shared.o -o $(STAR)

$(SRC)/mod_solvers.o: $(SRC)/mod_solvers.f90
	cd $(SRC); $(CC) -c mod_solvers.f90

$(SRC)/shooting_solution.o: $(SRC)/shooting_solution.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) -I$(SRC) -c shooting_solution.f90 mod_shared.o

$(SRC)/read_io.o: $(SRC)/read_io.f90 $(SRC)/mod_shared.o
	cd $(SRC); $(CC) -I$(SRC) -c read_io.f90 mod_shared.o

$(SRC)/mod_shared.o: $(SRC)/mod_shared.f90 
	cd $(SRC); $(CC) -c mod_shared.f90

clean: 
	rm src/*.o
	rm src/*.mod
