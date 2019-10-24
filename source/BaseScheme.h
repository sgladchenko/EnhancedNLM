#ifndef __BASESCHEME
#define __BASESCHEME

#include <mpi.h>

#define MPIComplex MPI_CXX_DOUBLE_COMPLEX
#define MPIWorld   MPI_COMM_WORLD

class BaseScheme
// The abstract, which determines the behaviour of each scheme:
// it must have methods of the initialisation from the binaries
// and it also must have the method of starting data processing
{
	public:
		// Initialise from the binaries
		virtual void Initialise() = 0;

		// Start data processing
		virtual void Process() = 0;

	protected:
		// The most significant info about the node
		int MyRank;
		int CommSize;

		// Ranks of the reighbour-nodes of this node
		int LeftNeighbour, RightNeighbour;
};

#endif