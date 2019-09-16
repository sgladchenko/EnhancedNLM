#include "inhomogeneities.h"
#include "constants.h"

#include <random>
#include <iostream>
#include <time.h>

using ph::Pi;
using ph::I;

// This function will be called only on master-node; after this call
// master should get (via broadcast, for example) all nodes array of harmonics

void generate_random_harmonics(int num_harmonics, double* harmonics)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, eps_mu / sqrt(num_harmonics));

	std::cout << "Generating harmonics:" << std::endl;

	generator.seed(time(NULL));

	for (int k = 1; k <= num_harmonics; ++k)
	{
		harmonics[k-1] = distribution(generator);
		std::cout << "k=" << k << ": " << harmonics[k-1] << std::endl;
	}

	std::cout << std::endl;
}

double mu(double X, double* harmonics, int num_harmonics)
{
	double tmp = 0;

	for (int k = 1; k <= num_harmonics; ++k)
	{
		tmp += harmonics[k-1] * sin(2.0*Pi*k*(X - Xmin)/L);
	}

	return ph::GF*sqrt(2.0)*n0 * (1.0 + tmp) * (1.0 - c2omega)*cubed(ph::hbar)*cubed(ph::c);
}