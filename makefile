CC=g++
CFLAGS=-o
LIBS = -fopenmp

ParallelIntegral: uebung1.cpp
	$(CC) $(CFLAGS) ParallelIntegral uebung1.cpp $(LIBS)
	
run: 
	./ParallelIntegral
