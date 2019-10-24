#include "Constants.h"
#include "InitialSpectra.h"
#include "Inhomogeneities.h"
#include "BaseBin.h"

#include <fstream>
#include <iostream>
using ph::Pi;

// Function to simplify writing of the binaries
template<typename T>
void bin_write(std::ofstream& stream, T* buf, int count)
{
	stream.write(reinterpret_cast<char*>(buf), count*sizeof(T));
}

class Generate
// The main class of generating the initial data for the calculational scheme
{
	public:
		void MakeRec();
		void MakeHarmonics();
		void MakeZ_init();

	private:
		// Functions to generate the initial conditions of the task
		double ee(double X, double E, InitialSpectra& Sp);
		double antiee(double X, double E, InitialSpectra& Sp);
		double ex(double X, double E, InitialSpectra& Sp);
		double antiex(double X, double E, InitialSpectra& Sp);
};

double Generate::ee(double X, double E, InitialSpectra& Sp)
{
	if (N_E > 1)
	{
		return Sp(E, 0) / (Sp(E, 0) + Sp(E, 1));
	}
	else
	{
		return ee0;
	}
}

double Generate::antiee(double X, double E, InitialSpectra& Sp)
{
	if (N_E > 1)
	{
		return Sp(E, 2) / (Sp(E, 2) + Sp(E, 3));
	}
	else
	{
		return antiee0;
	}
}

double Generate::ex(double X, double E, InitialSpectra& Sp)
{
	// To turn this off of we should merely choose eps_ex = 0
	return eps_ex * sin(2.0 * Pi * K * (X - Xmin) / (Xmax - Xmin));
}

double Generate::antiex(double X, double E, InitialSpectra& Sp)
{
	// To turn this off of we should merely choose eps_ex = 0
	return eps_ex * sin(2.0 * Pi * K * (X - Xmin) / (Xmax - Xmin));
}

void Generate::MakeRec()
{
	// Generate rec.bin file
	std::ofstream Rec(REC_BIN, std::ofstream::binary);

	// Auxiliary variables
	Complex ee, ex, xe, xx;
	InitialSpectra Sp;

	for (int x = 0; x < N_X; ++x)
	{
		// Current point of X-grid
		double X = Xmin + x*dX;

		for (int e = 0; e < N_E; ++e)
		{
			// Current point of energy in the spectra
			double E = Emin + e*dE;

			// Neutrino flow
			ee = this->ee(X, E, Sp);
			ex = this->ex(X, E, Sp);
			xe = std::conj(ex);
			xx = 1.0 - ee;
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				bin_write(Rec, &ee, 1);
				bin_write(Rec, &ex, 1);
				bin_write(Rec, &xe, 1);
				bin_write(Rec, &xx, 1);
			}

			// Antineutrino flow
			ee = this->antiee(X, E, Sp);
			ex = this->antiex(X, E, Sp);
			xe = std::conj(ex);
			xx = 1.0 - ee;
			for (int dz = 0; dz < 2; ++dz)
			{
				bin_write(Rec, &ee, 1);
				bin_write(Rec, &ex, 1);
				bin_write(Rec, &xe, 1);
				bin_write(Rec, &xx, 1);
			}
		}
	}

	Rec.close();	
}

void Generate::MakeZ_init()
{
	// Generate z.txt
	std::ofstream z(Z_INIT);
	z << 0;
	z.close();
}

void Generate::MakeHarmonics()
{
	// Generate harmonics.bin
	std::ofstream Harm(HARMONICS_BIN);

	double* Harmonics = new double[NUM_HARMONICS];
	GenerateRandomHarmonics(Harmonics);
	bin_write(Harm, Harmonics, NUM_HARMONICS);

	Harm.close();

	// Show these values
	std::cout << "Generated harmonics:" << std::endl;
	for (int k = 1; k <= NUM_HARMONICS; ++k)
	{
		std::cout << k << ": " << Harmonics[k-1] << std::endl;
	}
}

int main()
{
	Generate Gen;

	Gen.MakeRec();
	Gen.MakeZ_init();
	Gen.MakeHarmonics();
	
	return 0;
}