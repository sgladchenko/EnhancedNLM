#include "StandardScheme.h"
#include "AdiabaticScheme.h"

int main(int argc, char** argv)
{
	int CommSize, MyRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &CommSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

	// Scheme object
	AdiabaticScheme Scheme(MyRank, CommSize);

	// Initialise the grids and other data
	Scheme.Initialise();

	// Start data processing
	Scheme.Process();

	MPI_Finalize();
	return 0;
}