#include "matrix.h"
#include "pmatrix.h"
#include "InitialSpectra.h"
#include "constants.h"
#include "inhomogeneities.h"
#include "node.h"
#include <iostream>
#include <cmath>

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

void interchange_Ks(Matrix*  Kbuf,
					Complex* leftsend,
					Complex* rightsend,
					Complex* leftrecv,
					Complex* rightrecv,
					int left_neighbour,
					int right_neighbour,
					int Korder,
					int mycountX)
{
	// Some stuff to be used in the interchanging processes
	MPI_Request left_req, right_req;
	MPI_Status  left_st,  right_st;

	// Packing data, which's to be sent in futher

	for (int s = 0; s < N_E; ++s)
	{
		int x = 1;
		Kbuf[4*x*N_E + 4*s + 0].pack(leftsend + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].pack(leftsend + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].pack(leftsend + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].pack(leftsend + 16*N_E*(Korder-1) + 16*s + 3*4);

		x = mycountX;
		Kbuf[4*x*N_E + 4*s + 0].pack(rightsend + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].pack(rightsend + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].pack(rightsend + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].pack(rightsend + 16*N_E*(Korder-1) + 16*s + 3*4);

		//std::cout << leftsend[16*N_E*(Korder-1) + 16*s + 0*4] << std::endl;
	}

	// Interchanging the first coefficients of RK-cycle

	MPI_Isend(leftsend  + 16*N_E*(Korder-1), 16*N_E, MPIComplex, left_neighbour,  Korder, MPIWorld,  &left_req);
	MPI_Isend(rightsend + 16*N_E*(Korder-1), 16*N_E, MPIComplex, right_neighbour, Korder, MPIWorld, &right_req);

	MPI_Recv(leftrecv  + 16*N_E*(Korder-1),  16*N_E, MPIComplex, left_neighbour,  Korder, MPIWorld,  &left_st);
	MPI_Recv(rightrecv + 16*N_E*(Korder-1),  16*N_E, MPIComplex, right_neighbour, Korder, MPIWorld, &right_st);

	// Unpacking and saving at the grids

	for (int s = 0; s < N_E; ++s)
	{
		int x = 0;
		Kbuf[4*x*N_E + 4*s + 0].unpack(leftrecv + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].unpack(leftrecv + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].unpack(leftrecv + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].unpack(leftrecv + 16*N_E*(Korder-1) + 16*s + 3*4);

		x = mycountX + 1;
		Kbuf[4*x*N_E + 4*s + 0].unpack(rightrecv + 16*N_E*(Korder-1) + 16*s + 0*4);
		Kbuf[4*x*N_E + 4*s + 1].unpack(rightrecv + 16*N_E*(Korder-1) + 16*s + 1*4);
		Kbuf[4*x*N_E + 4*s + 2].unpack(rightrecv + 16*N_E*(Korder-1) + 16*s + 2*4);
		Kbuf[4*x*N_E + 4*s + 3].unpack(rightrecv + 16*N_E*(Korder-1) + 16*s + 3*4);	
	}
}

void interchange_PH_p(PMatrix* PH_p,
					  Complex* send_PH_p_buf,
					  Complex* recv_PH_p_buf,
					  int left_neighbour,
					  int right_neighbour,
					  int mycountX)
{
	// Packing all the interesting elements in the corresponding buffers

	for (int s = 0; s < N_E; ++s)
	{
		// Pack the rightest element

		int x = mycountX;
		PH_p[x*N_E + s].pack(send_PH_p_buf + s*4);
	}

	// Interchange of the values; send to the right neighbour and receive from the left one

	MPI_Request req; MPI_Status  stat;

	MPI_Isend(send_PH_p_buf, 4*N_E, MPIComplex, right_neighbour, 5, MPIWorld, &req);
	MPI_Recv(recv_PH_p_buf, 4*N_E, MPIComplex, left_neighbour, 5, MPIWorld, &stat);

	// Unpacking back the values from the receive buffer

	for (int s = 0; s < N_E; ++s)
	{
		// Unpack to the left element of the array

		int x = 0;
		PH_p[x*N_E + s].unpack(recv_PH_p_buf + s*4);
	}
}


void evaluate_adiabaticity(PMatrix* PH_p_old,
						   PMatrix* PH_p,
						   Complex* adiabaticity,
						   int mycountX)
{
	for (int x = 1; x <= mycountX; ++x)
	{
		for (int s = 0; s < N_E; ++s)
		{
			// The eigenvalues of the Hamiltonian here at this point (x, e)
			Complex eps_1, eps_2;

			// Normalise the Hamiltonian
			PH_p[x*N_E + s].normalise();

			// Obtain them from the PMatrix object
			PH_p[x*N_E + s].eigs(eps_1, eps_2);

			// Evaluate the derivative \nabla h
			PMatrix derivative;

			derivative = comega * (PH_p[x*N_E + s] - PH_p_old[x*N_E + s]) / dZ + 
						 somega * (PH_p[x*N_E + s] - PH_p[(x-1)*N_E + s]) / dX;

			// Make the skew product of PH_p and derivate
			PMatrix skew;

			skew = PH_p[x*N_E + s].skew(derivative);

			// Evaluate the factor and save it on the grid

			adiabaticity[(x-1)*N_E + s] = std::abs(eps_1 - eps_2) / skew.norm();
			//adiabaticity[(x-1)*N_E + s] = 5;
		}
	}
}

// RK-coefficients

void evaluate_K1(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 PMatrix* PH_p,
				 Matrix* K1) // the output buffer
{
	for (int x = 1; x <= mycountX; ++x)
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
		H_nn_p *= dE * mu(X + tanomega * Z, harmonics, NUM_HARMONICS);
		H_nn_m *= dE * mu(X - tanomega * Z, harmonics, NUM_HARMONICS);

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

			if (AD)
			{
				PH_p[x*N_E + s].from_matrix(H_p);
			}

			// Evaluate K1

			K1[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s]);     // p
			K1[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s]);     // m
			K1[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s]); // antip
			K1[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s]); // antim
		}
	}
}

void evaluate_K2(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K1,
				 Matrix* K2) // the output buffer
{
	for (int x = 1; x <= mycountX; ++x)
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
		H_nn_p *= dE * mu(X + tanomega * (Z + 0.5*dZ), harmonics, NUM_HARMONICS);
		H_nn_m *= dE * mu(X - tanomega * (Z + 0.5*dZ), harmonics, NUM_HARMONICS);

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

void evaluate_K3(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K2,
				 Matrix* K3) // the output buffer
{
	for (int x = 1; x <= mycountX; ++x)
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
		H_nn_p *= dE * mu(X + tanomega * (Z + 0.5*dZ), harmonics, NUM_HARMONICS);
		H_nn_m *= dE * mu(X - tanomega * (Z + 0.5*dZ), harmonics, NUM_HARMONICS);

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

			K3[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 0]*dZ);     // p
			K3[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 1]*dZ);     // m
			K3[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 2]*dZ); // antip
			K3[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + 0.5*K2[4*x*N_E + 4*s + 3]*dZ); // antim
		}
	}
}

void evaluate_K4(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K3,
				 Matrix* K4) // the output buffer
{
	for (int x = 1; x <= mycountX; ++x)
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
		H_nn_p *= dE * mu(X + tanomega * (Z + dZ), harmonics, NUM_HARMONICS);
		H_nn_m *= dE * mu(X - tanomega * (Z + dZ), harmonics, NUM_HARMONICS);

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

			K4[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + K3[4*x*N_E + 4*s + 0]*dZ);     // p
			K4[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + K3[4*x*N_E + 4*s + 1]*dZ);     // m
			K4[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + K3[4*x*N_E + 4*s + 2]*dZ); // antip
			K4[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + K3[4*x*N_E + 4*s + 3]*dZ); // antim
		}
	}
}
