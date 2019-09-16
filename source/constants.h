#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <math.h>
#include <complex>
typedef std::complex<double> Complex;

#define cubed(x) ((x)*(x)*(x))
#define squared(x) ((x)*(x))

/* PHYSICAL CONSTANTS */

namespace ph
{
	const double Pi = 3.1415926535897932384626433832795;
  	const double deg = Pi / 180;
  	const Complex I(0,1); 

	const double cm   = 1.0;
	const double sec  = 1.0;
	const double gram = 1.0;
	const double erg  = 1.0;

	const double eV   = 1.602e-12 * erg;
	const double MeV  = eV * 1e6;

	const double GF   = 1.166e-11 / (MeV*MeV); // Fermi constant
	const double hbar = 1.0546e-27 * erg * sec;
	const double c    = 2.9979e+10 * cm / sec;

	const double km = 1e5 * cm;
}

/* INHOMOGENEITY */

const int    NUM_HARMONICS = 25;
const double eps_ex = 0.0; // amplitude of inhomogenity in initial conditions
const double eps_mu = 0.1; // inhomogenetity in luminosity

/* MAIN CONSTANTS */

const double alpha   = 1;                    // the ratio of nu/antinu concentrations
const double ee0     = 0.666;                // values for spectra
const double xx0     = 1 - ee0;
const double antiee0 = 0.666;
const double antixx0 = 1 - antiee0;

const double eta   = 1.0;                    // the hierarchy of masses (either +1 or -1)
const double omega = 30.0 * ph::deg;         // the angle of both beams
const double theta = 9.0  * ph::deg;        // the lovely mixing angle
const double dm2   = 2.5e-3 * ph::eV*ph::eV; // the diff. between squared masses
const double n0    = 5.0e28 / cubed(ph::cm); // the constant level of the concentration, the 'luminosity'

const double c2omega  = cos(2*omega);
const double comega   = cos(omega);
const double tanomega = sin(omega) / cos(omega);
const double cotomega = cos(omega) / sin(omega);
const double c2theta  = cos(2*theta);
const double s2theta  = sin(2*theta);

/* E-MESH */

const double Emin = 1.0  * ph::MeV;      // reducing troubles which may occur with divisions by zero etc.
const double Emax = 50.0 * ph::MeV;      // the limit of 'infinite' integral
const int    N_E  = 10;                   // count of segments in energy spectra
const double dE   = (Emax - Emin) / N_E; // the step of energy mesh

/* X-MESH */

const double Xmin = 0.0  * ph::km;
const double Xmax = 50.0 * ph::km;
const double L    = Xmax - Xmin; 
const int    N_X  = 200;
const double dX   = Xmax / N_X;

/* Z-MESH OF THIS PART OF CALCULATIONS */

// the distance which will be observed in these calculations
// so if the initial point is Z0, the finite one will be Z0 + Z_displacement

const double Z_displacement = 50 * ph::km;

const int N_Z     = 250;  // the points which will be saved in files
const int STEP_Z  = 10; // the period of the mesh -- that means every STEP_Z'th point of real mesh will be placed in the mesh in output
const int COUNT_Z = N_Z * STEP_Z;  // the full number of points in real mesh, without the initial points
const double dZ   = Z_displacement / COUNT_Z;

#define NORM true

#endif