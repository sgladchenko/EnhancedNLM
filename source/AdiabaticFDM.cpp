#include "AdiabaticScheme.h"

void AdiabaticFDM(Layer* Rho_prev, Layer* Rho_here, Layer* Rho_next_0, Layer* Rho_next_1, Unit* H_0,
				  double* Spec_nu,
				  double* Spec_antinu,
				  double* Harmonics,
				  int     MyCountX,
				  double  Z,
				  double  Xleft);
{
	
}

// My army of the auxiliary functions, used in the FDM-scheme above

Complex NormDifferenceLayers(Layer* layer1, Layer* layer2, int MyCountX) // of course, all the norms are squared
{
	Complex Norm = 0;

	for (int x = 0; x < MyCountX; ++x)
	{
		Vector dv = layer1(x) - layer2(x);

		Norm += dv.norm();
	}

	return Norm;
}
