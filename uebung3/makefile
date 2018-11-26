LOCDIR = /usr/include/
DIRS = petsc slepc mpi/
LINKERFLAG = -lpetsc -lslepc -lm -lmpi
CC = g++

UEBUNG3_PC_POOL:
	${CC} -o uebung3 -I${LOCDIR}petsc -I${LOCDIR}slepc -I${LOCDIR}mpi/ ${LINKERFLAG} lambda.cpp hamilton.cpp diagonalize.cpp main.cpp
