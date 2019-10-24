#include "Matrix.h"
#include "InitialSpectra.h"
#include "Constants.h"
#include "Inhomogeneities.h"

#include <cmath>
#include <mpi.h>
#include <iostream>

// These functions are to implement the evaluation of the K1..K4
// coefficients of the RK-cycle

// We're adding them due to we have some new purposes as to 
// analise decoherence efffects, to examine what's happening with
// observed instabilities in the adiabatic limit.

// It would be a matter of great efforts, if we decided to apply new
// updates in the present code; hence the best decision is to realise
// this new code apart from the old one, and afterwards employ this in further,
// involving links to this module.

// imaginary unit and Pi
using ph::I;
using ph::Pi;

#define MPIComplex MPI_CXX_DOUBLE_COMPLEX
#define MPIWorld   MPI_COMM_WORLD

void Interchange_Ks(Matrix*  Kbuf,
					Complex* LeftSend_buf,
					Complex* RightSend_buf,
					Complex* LeftRecv_buf,
					Complex* RightRecv_buf,
					int      LeftNeighbour,
					int      RightNeighbour,
					int      Korder,
					int      MyCountX)
{
	// Some stuff to be used in the interchanging processes
	MPI_Request left_req, right_req;
	MPI_Status  left_st,  right_st;

	// Packing data, which's to be sent in futher

	for (int s = 0; s < N_E; ++s)
	{
		int x = 1;
		Kbuf[4*x*N_E + 4*s + 0].pack(LeftSend_buf + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].pack(LeftSend_buf + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].pack(LeftSend_buf + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].pack(LeftSend_buf + 16*N_E*(Korder-1) + 16*s + 3*4);

		x = MyCountX;
		Kbuf[4*x*N_E + 4*s + 0].pack(RightSend_buf + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].pack(RightSend_buf + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].pack(RightSend_buf + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].pack(RightSend_buf + 16*N_E*(Korder-1) + 16*s + 3*4);

		//std::cout << leftsend[16*N_E*(Korder-1) + 16*s + 0*4] << std::endl;
	}

	// Interchanging the first coefficients of RK-cycle

	MPI_Isend(LeftSend_buf  + 16*N_E*(Korder-1), 16*N_E, MPIComplex, LeftNeighbour,  Korder, MPIWorld,  &left_req);
	MPI_Isend(RightSend_buf + 16*N_E*(Korder-1), 16*N_E, MPIComplex, RightNeighbour, Korder, MPIWorld, &right_req);

	MPI_Recv(LeftRecv_buf  + 16*N_E*(Korder-1),  16*N_E, MPIComplex, LeftNeighbour,  Korder, MPIWorld,  &left_st);
	MPI_Recv(RightRecv_buf + 16*N_E*(Korder-1),  16*N_E, MPIComplex, RightNeighbour, Korder, MPIWorld, &right_st);

	// Unpacking and saving at the grids

	for (int s = 0; s < N_E; ++s)
	{
		int x = 0;
		Kbuf[4*x*N_E + 4*s + 0].unpack(LeftRecv_buf + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].unpack(LeftRecv_buf + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].unpack(LeftRecv_buf + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].unpack(LeftRecv_buf + 16*N_E*(Korder-1) + 16*s + 3*4);

		x = MyCountX + 1;
		Kbuf[4*x*N_E + 4*s + 0].unpack(RightRecv_buf + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].unpack(RightRecv_buf + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].unpack(RightRecv_buf + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].unpack(RightRecv_buf + 16*N_E*(Korder-1) + 16*s + 3*4);	
	}
}

// RK-coefficients

void Evaluate_K1(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K1) // the output buffer
{
	for (int x = 1; x <= MyCountX; ++x)
	{
		// Current X-coordinate

		double X = Xleft + (x - 1)*dX;

		// Integrals for Hamiltonians of zeta = +/-1 and for this X, from the previous iteration

		Matrix H_nn_p;
		Matrix H_nn_m;

		// integrating...
		for (int e = 0; e < N_E; ++e)
		{
			H_nn_p += Spec_nu[e]*old_m[x*N_E + e] - Spec_antinu[e]*old_antim[x*N_E + e];
			H_nn_m += Spec_nu[e]*old_p[x*N_E + e] - Spec_antinu[e]*old_antip[x*N_E + e];
		}

		// the collective parts of Hamiltonians are
		H_nn_p *= dE * mu(X + tanomega * Z, Harmonics);
		H_nn_m *= dE * mu(X - tanomega * Z, Harmonics);

		for (int s = 0; s < N_E; ++s)
		{
			// FDM-expressions for X-derivatives from the previous iteration
			Matrix Dp; Matrix Dantip;
			Matrix Dm; Matrix Dantim;

			Dp     = (0.5 / dX) * (old_p[(x+1)*N_E + s]     - old_p[(x-1)*N_E + s]);
			Dm     = (0.5 / dX) * (old_m[(x+1)*N_E + s]     - old_m[(x-1)*N_E + s]);
			Dantip = (0.5 / dX) * (old_antip[(x+1)*N_E + s] - old_antip[(x-1)*N_E + s]);
			Dantim = (0.5 / dX) * (old_antim[(x+1)*N_E + s] - old_antim[(x-1)*N_E + s]);

			// Hamiltonians, at these X and E
			Matrix H_p(H_0_mesh[s] + H_nn_p);
			Matrix H_m(H_0_mesh[s] + H_nn_m);
			Matrix H_antip(-1.0*H_0_mesh[s] + H_nn_p);
			Matrix H_antim(-1.0*H_0_mesh[s] + H_nn_m);

			// Save the PMatrix-formed Hamiltonians in order to interchange them in the following
			// (only if we need to evaluate the adiabaticity factor)

			// Evaluate K1
			K1[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s]);     // p
			K1[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s]);     // m
			K1[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s]); // antip
			K1[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s]); // antim
		}
	}
}

void Evaluate_K2(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K1,
				 Matrix* K2) // the output buffer
{
	for (int x = 1; x <= MyCountX; ++x)
	{
		// Current X-coordinate

		double X = Xleft + (x - 1)*dX;

		// Integrals for Hamiltonians of dzeta = +/-1 and for this X, from the previous iteration

		Matrix H_nn_p;
		Matrix H_nn_m;

		// integrating...
		for (int e = 0; e < N_E; ++e)
		{
			H_nn_p += Spec_nu[e]*    (old_m[x*N_E + e]     + 0.5*K1[4*x*N_E + 4*e + 1]*dZ)
					- Spec_antinu[e]*(old_antim[x*N_E + e] + 0.5*K1[4*x*N_E + 4*e + 3]*dZ);

			H_nn_m += Spec_nu[e]*    (old_p[x*N_E + e]     + 0.5*K1[4*x*N_E + 4*e + 0]*dZ)
					- Spec_antinu[e]*(old_antip[x*N_E + e] + 0.5*K1[4*x*N_E + 4*e + 2]*dZ);
		}

		// the collective parts of Hamiltonians are
		H_nn_p *= dE * mu(X + tanomega * (Z + 0.5*dZ), Harmonics);
		H_nn_m *= dE * mu(X - tanomega * (Z + 0.5*dZ), Harmonics);

		for (int s = 0; s < N_E; ++s)
		{
			// FDM-expressions for X-derivatives from the previous iteration
			Matrix Dp; Matrix Dantip;
			Matrix Dm; Matrix Dantim;

			Dp     = (0.5 / dX) * ((old_p[(x+1)*N_E + s] + 0.5*K1[4*(x+1)*N_E + 4*s + 0]*dZ)
								 - (old_p[(x-1)*N_E + s] + 0.5*K1[4*(x-1)*N_E + 4*s + 0]*dZ));

			Dm     = (0.5 / dX) * ((old_m[(x+1)*N_E + s] + 0.5*K1[4*(x+1)*N_E + 4*s + 1]*dZ)
								-  (old_m[(x-1)*N_E + s] + 0.5*K1[4*(x-1)*N_E + 4*s + 1]*dZ));

			Dantip = (0.5 / dX) * ((old_antip[(x+1)*N_E + s] + 0.5*K1[4*(x+1)*N_E + 4*s + 2]*dZ)
								-  (old_antip[(x-1)*N_E + s] + 0.5*K1[4*(x-1)*N_E + 4*s + 2]*dZ));

			Dantim = (0.5 / dX) * ((old_antim[(x+1)*N_E + s] + 0.5*K1[4*(x+1)*N_E + 4*s + 3]*dZ)
								-  (old_antim[(x-1)*N_E + s] + 0.5*K1[4*(x-1)*N_E + 4*s + 3]*dZ));

			// Hamiltonians, at these X and E
			Matrix H_p(H_0_mesh[s] + H_nn_p);
			Matrix H_m(H_0_mesh[s] + H_nn_m);
			Matrix H_antip(-1.0*H_0_mesh[s] + H_nn_p);
			Matrix H_antim(-1.0*H_0_mesh[s] + H_nn_m);

			// Evaluate K2
			K2[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 0]*dZ);     // p
			K2[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 1]*dZ);     // m
			K2[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 2]*dZ); // antip
			K2[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 3]*dZ); // antim
		}
	}
}

void Evaluate_K3(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K2,
				 Matrix* K3) // the output buffer
{
	for (int x = 1; x <= MyCountX; ++x)
	{
		// Current X-coordinate

		double X = Xleft + (x - 1)*dX;

		// Integrals for Hamiltonians of dzeta = +/-1 and for this X, from the previous iteration

		Matrix H_nn_p;
		Matrix H_nn_m;

		// integrating...
		for (int e = 0; e < N_E; ++e)
		{
			H_nn_p += Spec_nu[e]*    (old_m[x*N_E + e]     + 0.5*K2[4*x*N_E + 4*e + 1]*dZ)
					- Spec_antinu[e]*(old_antim[x*N_E + e] + 0.5*K2[4*x*N_E + 4*e + 3]*dZ);

			H_nn_m += Spec_nu[e]*    (old_p[x*N_E + e]     + 0.5*K2[4*x*N_E + 4*e + 0]*dZ)
					- Spec_antinu[e]*(old_antip[x*N_E + e] + 0.5*K2[4*x*N_E + 4*e + 2]*dZ);
		}

		// the collective parts of Hamiltonians are
		H_nn_p *= dE * mu(X + tanomega * (Z + 0.5*dZ), Harmonics);
		H_nn_m *= dE * mu(X - tanomega * (Z + 0.5*dZ), Harmonics);

		for (int s = 0; s < N_E; ++s)
		{
			// FDM-expressions for X-derivatives from the previous iteration
			Matrix Dp; Matrix Dantip;
			Matrix Dm; Matrix Dantim;

			Dp     = (0.5 / dX) * ((old_p[(x+1)*N_E + s] + 0.5*K2[4*(x+1)*N_E + 4*s + 0]*dZ)
								 - (old_p[(x-1)*N_E + s] + 0.5*K2[4*(x-1)*N_E + 4*s + 0]*dZ));

			Dm     = (0.5 / dX) * ((old_m[(x+1)*N_E + s] + 0.5*K2[4*(x+1)*N_E + 4*s + 1]*dZ)
								-  (old_m[(x-1)*N_E + s] + 0.5*K2[4*(x-1)*N_E + 4*s + 1]*dZ));

			Dantip = (0.5 / dX) * ((old_antip[(x+1)*N_E + s] + 0.5*K2[4*(x+1)*N_E + 4*s + 2]*dZ)
								-  (old_antip[(x-1)*N_E + s] + 0.5*K2[4*(x-1)*N_E + 4*s + 2]*dZ));

			Dantim = (0.5 / dX) * ((old_antim[(x+1)*N_E + s] + 0.5*K2[4*(x+1)*N_E + 4*s + 3]*dZ)
								-  (old_antim[(x-1)*N_E + s] + 0.5*K2[4*(x-1)*N_E + 4*s + 3]*dZ));

			// Hamiltonians, at these X and E
			Matrix H_p(H_0_mesh[s] + H_nn_p);
			Matrix H_m(H_0_mesh[s] + H_nn_m);
			Matrix H_antip(-1.0*H_0_mesh[s] + H_nn_p);
			Matrix H_antim(-1.0*H_0_mesh[s] + H_nn_m);

			// Evaluate K3
			K3[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 0]*dZ);     // p
			K3[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 1]*dZ);     // m
			K3[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 2]*dZ); // antip
			K3[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 3]*dZ); // antim
		}
	}
}

void Evaluate_K4(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K3,
				 Matrix* K4) // the output buffer
{
	for (int x = 1; x <= MyCountX; ++x)
	{
		// Current X-coordinate

		double X = Xleft + (x - 1)*dX;

		// Integrals for Hamiltonians of dzeta = +/-1 and for this X, from the previous iteration

		Matrix H_nn_p;
		Matrix H_nn_m;

		// integrating...
		for (int e = 0; e < N_E; ++e)
		{
			H_nn_p += Spec_nu[e]*    (old_m[x*N_E + e]     + K3[4*x*N_E + 4*e + 1]*dZ)
					- Spec_antinu[e]*(old_antim[x*N_E + e] + K3[4*x*N_E + 4*e + 3]*dZ);

			H_nn_m += Spec_nu[e]*    (old_p[x*N_E + e]     + K3[4*x*N_E + 4*e + 0]*dZ)
					- Spec_antinu[e]*(old_antip[x*N_E + e] + K3[4*x*N_E + 4*e + 2]*dZ);
		}

		// the collective parts of Hamiltonians are
		H_nn_p *= dE * mu(X + tanomega * (Z + dZ), Harmonics);
		H_nn_m *= dE * mu(X - tanomega * (Z + dZ), Harmonics);

		for (int s = 0; s < N_E; ++s)
		{
			// FDM-expressions for X-derivatives from the previous iteration
			Matrix Dp; Matrix Dantip;
			Matrix Dm; Matrix Dantim;

			Dp     = (0.5 / dX) * ((old_p[(x+1)*N_E + s] + K3[4*(x+1)*N_E + 4*s + 0]*dZ)
								 - (old_p[(x-1)*N_E + s] + K3[4*(x-1)*N_E + 4*s + 0]*dZ));

			Dm     = (0.5 / dX) * ((old_m[(x+1)*N_E + s] + K3[4*(x+1)*N_E + 4*s + 1]*dZ)
								-  (old_m[(x-1)*N_E + s] + K3[4*(x-1)*N_E + 4*s + 1]*dZ));

			Dantip = (0.5 / dX) * ((old_antip[(x+1)*N_E + s] + K3[4*(x+1)*N_E + 4*s + 2]*dZ)
								-  (old_antip[(x-1)*N_E + s] + K3[4*(x-1)*N_E + 4*s + 2]*dZ));

			Dantim = (0.5 / dX) * ((old_antim[(x+1)*N_E + s] + K3[4*(x+1)*N_E + 4*s + 3]*dZ)
								-  (old_antim[(x-1)*N_E + s] + K3[4*(x-1)*N_E + 4*s + 3]*dZ));

			// Hamiltonians, at these X and E
			Matrix H_p(H_0_mesh[s] + H_nn_p);
			Matrix H_m(H_0_mesh[s] + H_nn_m);
			Matrix H_antip(-1.0*H_0_mesh[s] + H_nn_p);
			Matrix H_antim(-1.0*H_0_mesh[s] + H_nn_m);

			// Evalaute K4
			K4[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + K3[4*x*N_E + 4*s + 0]*dZ);     // p
			K4[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + K3[4*x*N_E + 4*s + 1]*dZ);     // m
			K4[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + K3[4*x*N_E + 4*s + 2]*dZ); // antip
			K4[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + K3[4*x*N_E + 4*s + 3]*dZ); // antim
		}
	}
}
