CC = g++
FC = g95
CFLAGS = -O3 -Wall
LIBRARY_FLAGS = -L./lib -L/home/daschaich/gsl/lib -lgsl -lgslcblas -lm
FORTRAN_FLAGS = -L./lib -L/home/daschaich/gsl/lib -lgsl -lgslcblas -lm -lgfortran

EXECUTABLE = Simulation
OBJECT_FILES = HashTable.o Lattice.o Simulation.o
SIMPLE_FILES = HashTable.o Lattice_simple.o Simulation.o
WEIGHTED_FILES = HashTable.o Lattice_weighted.o Simulation.o
FFT_FILES = HashTable.o twodfft.o Lattice_FFT.o Simulation_FFT.o
4D_FILES = HashTable.o Lattice_4d.o Simulation_4d.o

main: HashTable.o Lattice.o Simulation.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECT_FILES) $(LIBRARY_FLAGS)

fast: HashTable.o Lattice.o Simulation.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECT_FILES) $(LIBRARY_FLAGS)

simple: HashTable.o Lattice_simple.o Simulation.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SIMPLE_FILES) $(LIBRARY_FLAGS)

weighted: HashTable.o Lattice_weighted.o Simulation.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(WEIGHTED_FILES) $(LIBRARY_FLAGS)

FFT: HashTable.o twodfft.o Lattice_FFT.o Simulation_FFT.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(FFT_FILES) $(FORTRAN_FLAGS)

4d: HashTable.o Lattice_4d.o Simulation_4d.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(4D_FILES) $(LIBRARY_FLAGS)

clean:
	rm -f $(EXECUTABLE) *.o

Simulation.o: Simulation.cpp Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Simulation.cpp

Simulation_FFT.o: Simulation_FFT.cpp Lattice_FFT.cpp Lattice_FFT.hh
	$(CC) $(CFLAGS) -c Simulation_FFT.cpp

Simulation_4d.o: Simulation_4d.cpp Lattice_4d.cpp Lattice_4d.hh
	$(CC) $(CFLAGS) -c Simulation_4d.cpp

Lattice.o: Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Lattice.cpp

Lattice_simple.o: Lattice_simple.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Lattice_simple.cpp

Lattice_weighted.o: Lattice_weighted.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Lattice_weighted.cpp

Lattice_FFT.o: Lattice_FFT.cpp Lattice_FFT.hh
	$(CC) $(CFLAGS) -c Lattice_FFT.cpp

Lattice_4d.o: Lattice_4d.cpp Lattice_4d.hh
	$(CC) $(CFLAGS) -c Lattice_4d.cpp

twodfft.o: twodfft.f90
	$(FC) $(CFLAGS) -c twodfft.f90

HashTable.o: HashTable.hh HashTable.cpp
	$(CC) $(CFLAGS) -c HashTable.cpp
