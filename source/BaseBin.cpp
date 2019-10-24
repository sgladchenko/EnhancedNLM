#include <fstream>
#include <iomanip>

#include "Segmentation.h"
#include "BaseBin.h"

using ph::km;
using ph::MeV;

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

BaseBin::BaseBin(int given_CommSize)
{
	CommSize = given_CommSize;

	CountXs = new int[CommSize];
	Obtain_CountXs(CommSize, CountXs);

	GatherCounts = new int[CommSize];
	Obtain_GatherCounts(CommSize, GatherCounts);

	GatherDisplacements = new int[CommSize];
	Obtain_GatherDisplacements(CommSize, GatherDisplacements);	
}

BaseBin::~BaseBin()
{
	delete [] CountXs;
	delete [] GatherCounts;
	delete [] GatherDisplacements;
}

// Input methods

void BaseBin::InputScatterBuffer(Complex* ScatterBuffer)
{
	std::ifstream rec;
	rec.open(REC_BIN, std::ifstream::binary);

	int index, length;

	index = 0;
	for (int r = 0; r < CommSize; ++r)
	{
		bin_read(rec, ScatterBuffer + index + N_E*4*4, CountXs[r]*N_E*4*4);
		index += (CountXs[r]+2)*N_E*4*4;
	}

	// and adjacent ones, for calculating the first coefficient of RK-cycle

	index  = -1; // index ScatterBuffer; used for writing
	length = 0;  // length of rec-steam which is already processed; used as an index in file

	for (int r = 0; r < CommSize; ++r)
	{
		// reading the first element in chunck
		bin_seek<Complex>(rec, length*N_E*4*4);

		if (length == 0)
		{
			// we need to put it to the chunck+2 of the last node
			bin_read(rec, ScatterBuffer + (N_X + 2*CommSize - 1)*N_E*4*4, N_E*4*4);
		}
		else
		{
			// regular element
			bin_read(rec, ScatterBuffer + index*N_E*4*4, N_E*4*4);
		}

		// reading the last element in chunck
		length += CountXs[r] - 1;
		index  += CountXs[r] + 3;

		bin_seek<Complex>(rec, length*N_E*4*4);

		if (length == N_X - 1)
		{
			// we need to put it to the chunck+2 of the first node
			bin_read(rec, ScatterBuffer, N_E*4*4);
		}
		else
		{
			// regular element
			bin_read(rec, ScatterBuffer + index*N_E*4*4, N_E*4*4);
		}

		// for the next element
		length += 1;
		index  -= 1;
	}

	rec.close();	
}

void BaseBin::InputZ_init(double& Z_init)
{
	std::ifstream z;
	z.open(Z_INIT);
	z >> Z_init;
	z.close();	
}

void BaseBin::InputHarmonics(double* Harmonics)
{
	std::ifstream harms;
	harms.open(HARMONICS_BIN, std::ifstream::binary);
	bin_read(harms, Harmonics, NUM_HARMONICS);
	harms.close();
}

// Output methods

void BaseBin::OutputRec(Complex* GatherBuffer)
{
	std::ofstream rec;
	rec.open(REC_BIN, std::ofstream::binary);

	for (int r = 0; r < CommSize; ++r)
	{
		int countX = GatherCounts[r] / (16*N_E);
		int dsp    = GatherDisplacements[r];

		for (int x = 0; x < countX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				bin_write(rec, GatherBuffer + dsp + 16*x*N_E + 16*e, 16);
			}
		}
	}

	rec.close();
}

void BaseBin::OutputLine(Complex* GatherBuffer)
{
	std::ofstream p, m, antip, antim;

	p.open(P_BIN, std::ios::app | std::ios::binary);
	m.open(M_BIN, std::ios::app | std::ios::binary);

	antip.open(ANTIP_BIN, std::ios::app | std::ios::binary);
	antim.open(ANTIM_BIN, std::ios::app | std::ios::binary);
	
	for (int r = 0; r < CommSize; ++r)
	{
		int countX = GatherCounts[r] / (16*N_E);
		int dsp    = GatherDisplacements[r];

		for (int x = 0; x < countX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				bin_write(p,     GatherBuffer + dsp + 16*x*N_E + 16*e + 0*4, 4);
				bin_write(m,     GatherBuffer + dsp + 16*x*N_E + 16*e + 1*4, 4);

				bin_write(antip, GatherBuffer + dsp + 16*x*N_E + 16*e + 2*4, 4);
				bin_write(antim, GatherBuffer + dsp + 16*x*N_E + 16*e + 3*4, 4);
			}
		}
	}

	p.close();
	m.close();
	antip.close();
	antim.close();
}

void BaseBin::OutputZ_final(double Z_final)
{
	std::ofstream z;
	z.open(Z_INIT);
	z << std::fixed;
	z << Z_final;
	z.close();
}

void BaseBin::MakeBins()
{
	std::ofstream p, m, antip, antim, ad;

	p.open(P_BIN, std::ios::binary);
	m.open(M_BIN, std::ios::binary);

	antip.open(ANTIP_BIN, std::ios::binary);
	antim.open(ANTIM_BIN, std::ios::binary);

	p.close();
	m.close();
	antip.close();
	antim.close();
}

void BaseBin::MakeZgrid()
{
	std::ofstream zgrid;
	zgrid.open(ZGRID, std::ios::binary);

	// Generate Z-grid file
	// Leave it empty!

	zgrid.close();
}

void BaseBin::MakeXgrid()
{
	std::ofstream xgrid;
	xgrid.open(XGRID, std::ios::binary);

	// Generate X-grid file

	for (int x = 0; x < N_X; ++x)
	{
		double X = (Xmin + x*dX) / ph::km;
		bin_write(xgrid, &X, 1);
	}

	xgrid.close();
}

void BaseBin::MakeEgrid()
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

void BaseBin::AddToZgrid(double Z_here)
{
	std::ofstream zgrid;
	zgrid.open(ZGRID, std::ios::app | std::ios::binary);

	// Add a new point to Z-grid

	double Z = Z_here / ph::km;
	bin_write(zgrid, &Z, 1);

	zgrid.close();
}
