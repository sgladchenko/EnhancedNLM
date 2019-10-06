#include "constants.h"
#include "segmentation.h"
#include "io.h"

#include <fstream>
#include <iomanip>

// Some templates to simplify writing/reading binaries

template<typename T>
void bin_read(std::ifstream& stream, T* buf, int count)
{
	stream.read(reinterpret_cast<char*>(buf), count*sizeof(T));
}

template<typename T>
void bin_seek(std::ifstream& stream, int count)
{
	stream.seekg(count*sizeof(T), std::ios::beg);
}

template<typename T>
void bin_write(std::ofstream& stream, T* buf, int count)
{
	stream.write(reinterpret_cast<char*>(buf), count*sizeof(T));
}

using ph::km;
using ph::MeV;

// Small extra function, used for initialization of the Z_init constant from another file

void Text_input_Z_init(double& Z_init)
{
	std::ifstream z;
	z.open(Z_FNAME);
	z >> Z_init;
	z.close();
}

// Small extra function, used for output Z_inal constant to file

void Text_output_Z_final(double Z_final)
{
	std::ofstream z;
	z.open(Z_FNAME);
	z << std::fixed;
	z << Z_final;
	z.close();
}

//  Function used for filling scatterv-buffer in master-node

void Bin_input_scatt_buffer(int size, Complex* scatt_buffer)
{
	std::ifstream rec;
	rec.open(REC_FNAME, std::ifstream::binary);

	int* countXs = new int[size];
	get_countXs(size, countXs);

	int index, length;

	// main chuncks belong to the nodes of ranks

	index = 0;
	for (int r = 0; r < size; ++r)
	{
		bin_read(rec, scatt_buffer + index + N_E*4*4, countXs[r]*N_E*4*4);
		index += (countXs[r]+2)*N_E*4*4;
	}

	// and adjacent ones, for calculating the first coefficient of RK-cycle

	index  = -1; // the index scatt_buffer; used for writing
	length = 0;  // the length of rec-steam which is already processed; used as an index in file

	for (int r = 0; r < size; ++r)
	{
		// reading the first element in chunck
		bin_seek<Complex>(rec, length*N_E*4*4);

		if (length == 0)
		{
			// we need to put it to the chunck+2 of the last node
			bin_read(rec, scatt_buffer + (N_X + 2*size - 1)*N_E*4*4, N_E*4*4);
		}
		else
		{
			// regular element
			bin_read(rec, scatt_buffer + index*N_E*4*4, N_E*4*4);
		}

		// reading the last element in chunck
		length += countXs[r] - 1;
		index  += countXs[r] + 3;

		bin_seek<Complex>(rec, length*N_E*4*4);

		if (length == N_X - 1)
		{
			// we need to put it to the chunck+2 of the first node
			bin_read(rec, scatt_buffer, N_E*4*4);
		}
		else
		{
			// regular element
			bin_read(rec, scatt_buffer + index*N_E*4*4, N_E*4*4);
		}

		// for the next element
		length += 1;
		index  -= 1;
	}

	delete [] countXs;

	rec.close();
}

// Function fills rec.bin file for following calculations

void Bin_output_rec(int size, int* counts, int* displacements, Complex* gath_buffer)
{
	std::ofstream rec;
	rec.open(REC_FNAME, std::ofstream::binary);

	for (int r = 0; r < size; ++r)
	{
		int countX = counts[r] / (16*N_E);
		int dsp    = displacements[r];

		for (int x = 0; x < countX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				bin_write(rec, gath_buffer + dsp + 16*x*N_E + 16*e, 16);
			}
		}
	}

	rec.close();
}

// Function fills harmonics array on master node

void Bin_input_harmonics(int num_harmonics, double* harmonics)
{
	std::ifstream harms;
	harms.open(HARMONICS_FILE, std::ifstream::binary);

	bin_read(harms, harmonics, num_harmonics);

	harms.close();
}

// Function used for binary output of one single line of fixed z variable

void make_bin_files()
{
	std::ofstream p, m, antip, antim, ad;

	p.open(P_OUTBIN, std::ios::binary);
	m.open(M_OUTBIN, std::ios::binary);

	antip.open(ANTIP_OUTBIN, std::ios::binary);
	antim.open(ANTIM_OUTBIN, std::ios::binary);

	ad.open(AD_OUT, std::ios::binary);

	p.close();
	m.close();
	antip.close();
	antim.close();

	ad.close();
}

void Bin_output_line(int size, int* counts, int* displacements, Complex* gath_buffer)
{
	std::ofstream p, m, antip, antim;

	p.open(P_OUTBIN, std::ios::app | std::ios::binary);
	m.open(M_OUTBIN, std::ios::app | std::ios::binary);

	antip.open(ANTIP_OUTBIN, std::ios::app | std::ios::binary);
	antim.open(ANTIM_OUTBIN, std::ios::app | std::ios::binary);
	
	for (int r = 0; r < size; ++r)
	{
		int countX = counts[r] / (16*N_E);
		int dsp    = displacements[r];

		for (int x = 0; x < countX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				bin_write(p,     gath_buffer + dsp + 16*x*N_E + 16*e + 0*4, 4);
				bin_write(m,     gath_buffer + dsp + 16*x*N_E + 16*e + 1*4, 4);

				bin_write(antip, gath_buffer + dsp + 16*x*N_E + 16*e + 2*4, 4);
				bin_write(antim, gath_buffer + dsp + 16*x*N_E + 16*e + 3*4, 4);
			}
		}
	}

	p.close();
	m.close();
	antip.close();
	antim.close();
}

// Updates of files of grids (they're contains arrays of points of each axis)

void make_gridx()
{
	std::ofstream xgrid;
	xgrid.open(XGRID, std::ios::binary);

	// Generate x-grid file

	for (int x = 0; x < N_X; ++x)
	{
		double X = (Xmin + x*dX) / ph::km;
		bin_write(xgrid, &X, 1);
	}

	xgrid.close();
}

void make_gridE()
{
	std::ofstream egrid;
	egrid.open(EGRID, std::ios::binary);

	// Generate E-grid file

	for (int e = 0; e < N_E; ++e)
	{
		double E = (Emin + e*dE) / ph::MeV;
		bin_write(egrid, &E, 1);
	}

	egrid.close();
}

void make_gridz()
{
	std::ofstream zgrid;
	zgrid.open(ZGRID, std::ios::binary);

	// Generate z-grid file

	// Leave it empty

	zgrid.close();
}

void add_to_gridz(double Z_here)
{
	std::ofstream zgrid;

	zgrid.open(ZGRID, std::ios::app | std::ios::binary);

	// Add new point to z-grid

	double Z = Z_here / ph::km;

	bin_write(zgrid, &Z, 1);

	zgrid.close();
}

// Binary output of the grid of the adiabaticity factor

void Bin_output_ad(int size, int* ad_counts, int* ad_displacements, Complex* adiabaticity_buffer)
{
	std::ofstream ad_output;

	ad_output.open(AD_OUT, std::ios::app | std::ios::binary);

	for (int r = 0; r < size; ++r)
	{
		int countX = ad_counts[r] / N_E;
		int dsp    = ad_displacements[r];

		for (int x = 0; x < countX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				bin_write(ad_output, adiabaticity_buffer + dsp + x*N_E + e, 1);
			}
		}
	}

	ad_output.close();
}