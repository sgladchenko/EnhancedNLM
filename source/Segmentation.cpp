#include "Segmentation.h"
#include "Constants.h"

void Obtain_SegmentX(int MyRank, int CommSize, int& MyLeftX, int& MyRightX) // [MyLeftX, MyRightX] in [0; N_X-1]
{
	// Function gives the chunck of X cooridinates

	if (MyRank >= N_X % CommSize)
	{
		MyLeftX = (N_X % CommSize) * (N_X / CommSize + 1) + (MyRank - N_X % CommSize) * (N_X / CommSize);
		MyRightX = MyLeftX + N_X / CommSize - 1;
	}
	else
	{
		MyLeftX = MyRank * (N_X / CommSize + 1);
		MyRightX = MyLeftX + N_X / CommSize;
	}

	// Note: MyLeftX and MyRightX are the indecies of segments, which belong to the unit of rank!
	// That means MyRightX'th segment also belongs to particular node!
}

void Obtain_CountXs(int CommSize, int* CountXs)
{
	// Returns the array containing the countX's for each node in communicator
	
	int tmpLeftX, tmpRightX;

	for (int r = 0; r < CommSize; ++r)
	{
		Obtain_SegmentX(r, CommSize, tmpLeftX, tmpRightX);
		CountXs[r] = tmpRightX - tmpLeftX + 1;
	}
}

void Obtain_ScatterCounts(int CommSize, int* ScatterCounts)
{
	// Array contains the counts for each node in the communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int tmpLeftX, tmpRightX;

	for (int r = 0; r < CommSize; ++r)
	{
		Obtain_SegmentX(r, CommSize, tmpLeftX, tmpRightX);
		ScatterCounts[r] = (tmpRightX - tmpLeftX + 1 + 2) * 16*N_E;
	}
}

void Obtain_ScatterDisplacements(int CommSize, int* ScatterDisplacements)
{
	// Array contains the displacements for each node in the communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int tmpLeftX, tmpRightX;
	int tmp_sum = 0;

	for (int r = 0; r < CommSize; ++r)
	{
		Obtain_SegmentX(r, CommSize, tmpLeftX, tmpRightX);
		ScatterDisplacements[r] = tmp_sum * 16*N_E;
		tmp_sum += tmpRightX - tmpLeftX + 1 + 2;
	}
}

void Obtain_GatherCounts(int CommSize, int* GatherCounts)
{
	// Array contains the counts for each node in the communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int tmpLeftX, tmpRightX;

	for (int r = 0; r < CommSize; ++r)
	{
		Obtain_SegmentX(r, CommSize, tmpLeftX, tmpRightX);
		GatherCounts[r] = (tmpRightX - tmpLeftX + 1) * 16*N_E;
	}
}

void Obtain_GatherDisplacements(int CommSize, int* GatherDisplacements)
{
	// Array contains the displacements for each node in the communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int tmpLeftX, tmpRightX;
	int tmp_sum = 0;

	for (int r = 0; r < CommSize; ++r)
	{
		Obtain_SegmentX(r, CommSize, tmpLeftX, tmpRightX);
		GatherDisplacements[r] = tmp_sum * 16*N_E;
		tmp_sum += tmpRightX - tmpLeftX + 1;
	}
}