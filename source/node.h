#ifndef __NODE_H
#define __NODE_H

#include <mpi.h>

#define MPIComplex MPI_CXX_DOUBLE_COMPLEX
#define MPIWorld   MPI_COMM_WORLD

void node(int rank, int size);

#endif