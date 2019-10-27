#include "AdiabaticScheme.h"

using ph::I;
using ph::Pi;

void AdiabaticFDM(Layer* Rho_prev, Layer* Rho_here, Layer* Rho_next_0, Layer* Rho_next_1,
				  Layer* H_prev,   Layer* H_here,   Layer* H_next_0,   Layer* H_next_1,
				  Unit* H_0,
				  double* Spec_nu,
				  double* Spec_antinu,
				  double* Harmonics,
				  int     MyCountX,
				  double  Z,
				  double  Xleft)
{
	for (int x = 1; x <= MyCountX; ++x)
	{
		Vector Rho_km1_n = (*Rho_here)(x-1);
		Vector Rho_kp1_n = (*Rho_here)(x+1);

		Vector Rho_k_n = (*Rho_here)(x);
		Vector H_k_n   = (*H_here)(x);

		Vector term1 = -tanomega * dZ/dX * (Rho_kp1_n - Rho_km1_n);
		Vector term2 = -2.0*I*dZ/(ph::hbar*ph::c*comega) * commutator(H_k_n, Rho_k_n);

		term1.MultiplyByPrefactor();

		(*Rho_next_0)(x) = (*Rho_prev)(x) + term1 + term2;
	}
}

// My army of the auxiliary functions, used in the FDM-scheme above

double Indicator(int type_or_zeta)
{
	if (type_or_zeta == 0)
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}

int OppositeZeta(int zeta)
{
	if (zeta == 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void EvaluateHamiltonian(Layer* Rho_buf, Layer* H_buf, Unit* H_0,
			     		 double* Spec_nu,
				  		 double* Spec_antinu,
				  		 double* Harmonics,
				  		 int     MyCountX,
				  		 double  Z,
				  		 double  Xleft)
{
	for (int x = 0; x <= MyCountX+1; ++x)
	{
		Vector Rho_x = (*Rho_buf)(x); // Pick a particular Vector-object at this coordinate x
		Vector H_x;                   // The Hamiltonian's Vector-object at this coordinate x to be assembled

		double X = Xleft + (x - 1)*dX; // The coordinate

		for (int type = 0; type < 2; ++type)
		{
			double indicator_type = Indicator(type);

			for (int zeta = 0; zeta < 2; ++zeta)
			{
				int opposite_zeta  = OppositeZeta(zeta);
				int indicator_zeta = Indicator(zeta);

				// Evaluating the interaction part of the Hamiltonian
				Unit interaction;

				for (int e = 0; e < N_E; ++e)
				{
					interaction += Spec_nu[e]*Rho_x(e, 0, opposite_zeta) - Spec_antinu[e]*Rho_x(e, 1, opposite_zeta);
				}

				// Either (X + tanomega*Z) for zeta=0 or (X - tanomega*Z) for zeta=1
				interaction *= dE * mu(X  + indicator_zeta*tanomega * Z, Harmonics);

				for (int e = 0; e < N_E; ++e)
				{
					H_x(e, type, zeta) = indicator_type*H_0[e] + interaction;
				}
			}
		}

		(*H_buf)(x) = H_x;
	}
}

Complex NormDifferenceLayers(Layer* layer1, Layer* layer2, int MyCountX) // of course, all the norms are squared
{
	Complex Norm = 0;

	for (int x = 1; x <= MyCountX; ++x)
	{
		Vector dv = (*layer1)(x) - (*layer2)(x);
		Norm += dv.norm();
	}

	return Norm;
}

