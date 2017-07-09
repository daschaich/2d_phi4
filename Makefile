CC = g++
CFLAGS = -O3 -Wall
LIBRARY_FLAGS = -lgsl -lgslcblas -lm

EXECUTABLE = 2d_phi4
OBJECT_FILES = HashTable.o Lattice.o Simulation.o

main: HashTable.o Lattice.o Simulation.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECT_FILES) $(LIBRARY_FLAGS)

clean:
	rm -f $(EXECUTABLE) *.o

Simulation.o: Simulation.cpp Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Simulation.cpp

Lattice.o: Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Lattice.cpp

HashTable.o: HashTable.hh HashTable.cpp
	$(CC) $(CFLAGS) -c HashTable.cpp
