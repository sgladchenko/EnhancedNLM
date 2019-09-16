#include <fstream>
#include <iomanip>

#include "constants.h"
#include "io.h"
#include "InitialSpectra.h"

#include "inhomogeneities.h"

using ph::Pi;

template<typename T>
void bin_write(std::ofstream& stream, T* buf, int count)
{
	stream.write(reinterpret_cast<char*>(buf), count*sizeof(T));
}

double rho0_ee(double X, double E, InitialSpectra& Sp)
{
	return Sp(E, 0) / (Sp(E, 0) + Sp(E, 1));
	//return ee0;
}

double rho0_ex(double X, double E, InitialSpectra& Sp)
{
	return 0;
	//return eps_ex * sin(2.0 * Pi * K * (X - Xmin) / (Xmax - Xmin));
}


double rho0_antiee(double X, double E, InitialSpectra& Sp)
{
	return Sp(E, 2) / (Sp(E, 2) + Sp(E, 3));
	//return antiee0;
}

double rho0_antiex(double X, double E, InitialSpectra& Sp)
{
	return 0;
}

int main()
{
	// Generate rec.bin file

	std::ofstream rec(REC_FNAME, std::ofstream::binary);

	Complex ee, ex, xe, xx;
	InitialSpectra Sp;

	for (int x = 0; x < N_X; ++x)
	{
		// Current point of X-mesh

		double X = Xmin + x*dX;

		for (int e = 0; e < N_E; ++e)
		{
			// Current point of energy in spectra

			double E = Emin + e*dE;

			ee = rho0_ee(X, E, Sp);
			ex = rho0_ex(X, E, Sp);
			xe = std::conj(ex);
			xx = 1.0 - ee;

			for (int dz = 0; dz < 2; ++dz)
			{
				bin_write(rec, &ee, 1);
				bin_write(rec, &ex, 1);
				bin_write(rec, &xe, 1);
				bin_write(rec, &xx, 1);
			}

			ee = rho0_antiee(X, E, Sp);
			ex = rho0_antiex(X, E, Sp);
			xe = std::conj(ex);
			xx = 1.0 - ee;

			for (int dz = 0; dz < 2; ++dz)
			{
				bin_write(rec, &ee, 1);
				bin_write(rec, &ex, 1);
				bin_write(rec, &xe, 1);
				bin_write(rec, &xx, 1);
			}
		}
	}

	rec.close();

	// Generate z.txt

	std::ofstream z(Z_FNAME);
	z << 0;
	z.close();

	// Generate harmonics.bin

	std::ofstream harms(HARMONICS_FILE);

	double* harmonics = new double[NUM_HARMONICS];

	generate_random_harmonics(NUM_HARMONICS, harmonics);
	bin_write(harms, harmonics, NUM_HARMONICS);

	harms.close();

	return 0;
}