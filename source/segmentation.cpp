#include "segmentation.h"
#include "constants.h"

// Function gives the chunck of X cooridinates
void get_segment_X(int rank, int size, int& left, int& right) // [left, right] in [0; N_X-1]
{
	if (rank >= N_X%size)
	{
		left  = (N_X%size)*(N_X/size+1) + (rank - N_X%size)*(N_X/size);
		right = left + N_X/size - 1;
	}
	else
	{
		left  = rank*(N_X/size+1);
		right = left + N_X/size;
	}

	// Note: left and right are the indecies of segments which belong to the unit of rank!
	// That means right'th segment also belongs to particular unit!
	// Note: and meshes will contain the left points of each segment; so each unit will
	// have only left points and the amount of points in the whole mesh is N_X (the last will be N_X-1)
}

void get_countXs(int size, int* countXs)
{
	// Array contains the countX's for each node in communicator
	
	int temp_leftX, temp_rightX;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		countXs[r] = temp_rightX - temp_leftX + 1;
	}
}

void get_scatter_counts(int size, int* counts)
{
	// Array contains the counts for each node in communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int temp_leftX, temp_rightX;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		counts[r] = (temp_rightX - temp_leftX + 1 + 2) * 16*N_E;
	}
}

void get_scatter_displacements(int size, int* displacements)
{
	// Array contains the displacements for each node in communicator
	// while scattering
	// (in sizes of MPIComplex)
	
	int temp_leftX, temp_rightX, temp_sum;
	temp_sum = 0;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		displacements[r] = temp_sum * 16 * N_E;
		temp_sum += temp_rightX - temp_leftX + 1 + 2;
	}
}

// These functions are used for gathering small amounts of data
// (it appeared in update when we added gathering while calculations are being processed)

void get_line_counts(int size, int* counts)
{
	// Array contains the counts for each node in communicator
	// while gathering ONE line of constant z
	// (in sizes of MPIComplex)
	
	int temp_leftX, temp_rightX;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		counts[r] = (temp_rightX - temp_leftX + 1) * 16 * N_E;
	}
}

void get_line_displacements(int size, int* displacements)
{
	// Array contains the displacements for each node in communicator
	// while gathering ONE line of constant z
	// (in sizes of MPIComplex)
	
	int temp_leftX, temp_rightX, temp_sum;
	temp_sum = 0;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		displacements[r] = temp_sum * 16 * N_E;
		temp_sum += temp_rightX - temp_leftX + 1;
	}
}

void get_ad_counts(int size, int* ad_counts)
{
	// Similar functionality, but for the gathering of adiabaticity factors

	int temp_leftX, temp_rightX;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		ad_counts[r] = (temp_rightX - temp_leftX + 1) * N_E;
	}

}

void get_ad_displacements(int size, int* ad_displacements)
{
	// Similar functionality, but for the gathering of adiabaticity factors

	int temp_leftX, temp_rightX, temp_sum;
	temp_sum = 0;

	for (int r = 0; r < size; ++r)
	{
		get_segment_X(r, size, temp_leftX, temp_rightX);
		ad_displacements[r] = temp_sum * N_E;
		temp_sum += temp_rightX - temp_leftX + 1;
	}
}