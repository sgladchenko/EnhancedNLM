#include "BaseScheme.h"
#include "BaseBin.h"

#include "InitialSpectra.h"
#include "Segmentation.h"
#include "Constants.h"
#include "Matrix.h"
#include "Containers.h"

void AdiabaticFDM(Layer* Rho_prev, Layer* Rho_here, Layer* Rho_next_0, Layer* Rho_next_1, Unit* H_0,
				  double* Spec_nu,
				  double* Spec_antinu,
				  double* Harmonics,
				  int     MyCountX,
				  double  Z,
				  double  Xleft);

// Remark: Rho_next_0 is a kth approximation and Rho_next_1 is a (k+1)th one

class AdiabaticScheme: public BaseScheme
// The implementation of the numerical scheme with
// included adiabaticity
{
	public:
		AdiabaticScheme(int given_MyRank, int given_CommSize);

		void Initialise();
		void Process();

		~AdiabaticScheme();

	protected:
		// Energy spectra to be used in the integrals etc.
		InitialSpectra Sp;

		// X-chunck of this node
		int MyCountX, MyLeftX, MyRightX;

		// The particular values of the X coordinates
		double Xleft, Xright;

		// Initial Z-coordinate
		double Z_init;

		// Buffer of the initial scattering
		Complex* ScatterBuffer;

		// Buffer of the gathering (for saving a new line of data)
		Complex* GatherBuffer;

		// Array of the harmonics of inhomogeneity 
		double* Harmonics;

		// Something for the scattering
		int* ScatterCounts;
		int* ScatterDisplacements;

		// Something for the scattering
		int* GatherCounts;
		int* GatherDisplacements;

		// Data to be used in the calculational process
		Layer* Rho_prev;
		Layer* Rho_here;
		Layer* Rho_next_0;
		Layer* Rho_next_1;

		Layer Rhos[4]; // Real objects of the aliases above

		// Spectra buffers; used in integrals
		double* Spec_e;
		double* Spec_x;

		double* Spec_antie;
		double* Spec_antix;

		double* Spec_nu;
		double* Spec_antinu;

		// Buffer which contains the constant terms in Hamiltonians,
		// depends only on energy
		Unit* H_0;

		// The very start time stamp of the calculations
		double Start_time;
};