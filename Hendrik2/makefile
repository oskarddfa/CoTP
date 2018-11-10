CC=g++
CFLAGS= -Wall
LIBS = -lpetsc -lslepc -lm -pthread -lmpi -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/home/hendrik/Documents/ROOT/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
INCLUDE =  -I/usr/include/slepc -I/usr/include/petsc -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -pthread -std=c++11 -m64 -I/home/hendrik/Documents/ROOT/include
OBJFILES = uebung2.o Potentiale.o ROOTplot.o

ROOTLIBS = -L/home/hendrik/Documents/ROOT/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOTCFLAGS = -pthread -std=c++11 -m64 -I/home/hendrik/Documents/ROOT/include

TARGET = Schroedinger



all: $(TARGET)

$(TARGET): $(OBJFILES)

	$(CC) -o  $(TARGET) $(CFLAGS) $(ROOTCFLAGS) $^ $(ROOTLIBS) $(LIBS) 
	make clean
	
uebung2.o: uebung2.cpp 

	$(CC) -c $(CFLAGS) $^ $(INCLUDE)

Potentiale.o: Potentiale.cpp

	$(CC) $(CFLAGS) -c $<
	
ROOTplot.o: ROOTplot.cpp

	$(CC) $(CFLAGS) $(ROOTCFLAGS) -c $< 
run:
	./$(TARGET)
	
clean: 
	$(RM) $(OBJFILES)
	


