#include "Inhomogeneities.h"
#include "Constants.h"

#include <random>
#include <time.h>

using ph::Pi;
using ph::I;

// This function will be called only on master-node; after this call
// master should get (via broadcast, for example) all nodes array of harmonics

void GenerateRandomHarmonics(double* Harmonics)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, eps_mu / sqrt(NUM_HARMONICS));

	generator.seed(time(NULL));

	for (int k = 1; k <= NUM_HARMONICS; ++k)
	{
		Harmonics[k-1] = distribution(generator);
	}
}

// The main function describining the inhomogeneity in the luminocity
// of the supernova

double mu(double X, double* Harmonics)
{
	double tmp = 0;

	for (int k = 1; k <= NUM_HARMONICS; ++k)
	{
		tmp += Harmonics[k-1] * sin(2.0*Pi*k*(X - Xmin)/L);
	}

	return ph::GF*sqrt(2.0)*n0 * (1.0 + tmp) * (1.0 - c2omega)*cubed(ph::hbar)*cubed(ph::c);
}