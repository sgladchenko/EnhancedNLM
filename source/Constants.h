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
const double eps_ex = 0.0;       // Amplitude of inhomogenity in initial conditions (inhomogeneity of the off-diag. element)
const double eps_mu = 0.1;       // Inhomogenetity in the luminosity
const int    K      = 0;         // Might be used if we were studying instabilities as Duan did

/* MAIN CONSTANTS */

// Ratio of nu/antinu concentrations
const double alpha   = 1;
// Values for spectra (will be used in the monochromatic configuration N_E=1)
const double ee0     = 0.666;
const double xx0     = 1 - ee0;
const double antiee0 = 0.666;
const double antixx0 = 1 - antiee0;

// Hierarchy of masses (either +1 or -1)
const double eta   = 1.0;
// Angle of both beams
const double omega = 30.0 * ph::deg;
// Lovely mixing angle
const double theta = 9.0  * ph::deg; 
// Difference between the squared masses
const double dm2   = 2.5e-3 * ph::eV*ph::eV;
// Constant level of the concentration, the 'luminosity'
const double n0    = 1.0e28 / cubed(ph::cm);

// Different functions of the angles, which take place in the equation
const double c2omega  = cos(2*omega);
const double comega   = cos(omega);
const double somega   = sin(omega);
const double tanomega = sin(omega) / cos(omega);
const double cotomega = cos(omega) / sin(omega);
const double c2theta  = cos(2*theta);
const double s2theta  = sin(2*theta);

/* E-GRID */

// This's not zero, resolving troubles which may occur with divisions by zero etc.
const double Emin = 1.0  * ph::MeV;
// Limit of the 'infinite' integral
const double Emax = 50.0 * ph::MeV;
// Count of the segments in the energy spectrum
const int    N_E  = 20;
// Step of the energy grid
const double dE   = (Emax - Emin) / N_E;

/* X-GRID */

// Left point of the x-grid (normally it's zero)
const double Xmin = 0.0  * ph::km;
// Right one
const double Xmax = 50.0 * ph::km;
// The length of the grid
const double L = Xmax - Xmin; 
// Number of points of this grid
const int N_X  = 300;
// Step of this grid
const double dX = Xmax / N_X;

/* Z-GRID OF THIS PART OF CALCULATIONS */

// The distance which will be processed in these calculations
// i.e. if the initial point is Z0, the finite one will be Z0 + Z_displacement
const double Z_displacement = 50 * ph::km;
// Number of points to be saved in the binaries
const int N_Z = 500;
// Period of the z-grid: i.e. every STEP_Z'th point of the real grid will be placed in the output
const int STEP_Z = 5;
// Full number of points in the real z-grid, without the initial ones
const int COUNT_Z = N_Z * STEP_Z;            
// Step of the grid 
const double dZ = Z_displacement / COUNT_Z;

/* SOME OTHER FLAGS */

// To normalise the density matrices after each iteration or not (recovering their hermitance and trace)
#define NORM true

#endif