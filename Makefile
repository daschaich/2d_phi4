CC = g++
CFLAGS = -O3 -Wall
LIBRARY_FLAGS = -lgsl -lgslcblas -lm

EXECUTABLE = 2d_phi4
OBJECT_FILES = HashTable.o Lattice.o 2d_phi4.o

main: HashTable.o Lattice.o 2d_phi4.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECT_FILES) $(LIBRARY_FLAGS)

clean:
	rm -f $(EXECUTABLE) *.o

2d_phi4.o: 2d_phi4.cpp Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c 2d_phi4.cpp

Lattice.o: Lattice.cpp Lattice.hh
	$(CC) $(CFLAGS) -c Lattice.cpp

HashTable.o: HashTable.hh HashTable.cpp
	$(CC) $(CFLAGS) -c HashTable.cpp
