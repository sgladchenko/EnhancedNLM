#include "node.h"
#include "segmentation.h"
#include "io.h"
#include "matrix.h"
#include "InitialSpectra.h"
#include "inhomogeneities.h"

#include <iostream>
#include <iomanip>

// Pauli matrices which will be used in Hamiltonians
Matrix tau1(0, 1, 1, 0);
Matrix tau3(1, 0, 0, -1);

// imaginary unit and Pi
using ph::I;
using ph::Pi;

void node(int rank, int size)
{
	// BASICS //

	// Spectra object

	InitialSpectra Sp;

	// Initialize my personal chunck of X-points

	int mycountX, myleftX, myrightX;
	get_segment_X(rank, size, myleftX, myrightX);
	mycountX  = myrightX - myleftX + 1;

	double Xleft  = Xmin + myleftX*dX;
	double Xright = Xmin + myrightX*dX;

	// Ranks of neighbours of this rank-node

	int left_neighbour, right_neighbour;

	if (rank == 0)
	{
		left_neighbour  = size - 1;
		right_neighbour = 1;
	}
	else if (rank == size - 1)
	{
		left_neighbour  = size - 2;
		right_neighbour = 0;
	}
	else
	{
		left_neighbour  = rank - 1;
		right_neighbour = rank + 1;
	}

	// INITIALIZATION: READING FILES //

	Complex* scatt_buffer;
	double Z_init;

	if (rank == 0)
	{
		scatt_buffer = new Complex[(N_X + 2*size) * 4 * 4 * N_E];
		Bin_input_scatt_buffer(size, scatt_buffer);
		Text_input_Z_init(Z_init);
	}
	else
	{
		scatt_buffer = new Complex[(mycountX + 2) * 4 * 4 * N_E];
	}

	// Local buffer used for gathering

	Complex* gath_buffer;

	if (rank == 0)
	{
		gath_buffer = new Complex[N_X * 16*N_E];
	}
	else
	{
		gath_buffer = new Complex[mycountX * 16*N_E];
	}

	// INITIALIZATION: SCATTERING //

	// Some structure for sending the entire spectra, with all elements for all dzeta and nu/antinu

	MPI_Datatype MPISpectra;
	MPI_Type_contiguous(4 * 4 * N_E, MPIComplex, &MPISpectra);
	MPI_Type_commit(&MPISpectra);

	int* counts        = new int[size];
	int* displacements = new int[size];

	get_scatter_counts(size, counts);
	get_scatter_displacements(size, displacements);

	MPI_Scatterv(scatt_buffer, counts, displacements, MPIComplex, scatt_buffer, (mycountX + 2) * 16*N_E, MPIComplex, 0, MPIWorld);
	MPI_Bcast(&Z_init, 1, MPI_DOUBLE, 0, MPIWorld);

	// INITIALIZATION: GENERATING RANDOM NOISE //

	double* harmonics;

	if (NUM_HARMONICS > 0)
	{
		harmonics = new double[NUM_HARMONICS];

		if (rank == 0)
		{
			Bin_input_harmonics(NUM_HARMONICS, harmonics);

			for (int k = 1; k <= NUM_HARMONICS; ++k)
			{
				std::cout << harmonics[k-1] << std::endl;
			}
		}

		MPI_Bcast(harmonics, NUM_HARMONICS, MPI_DOUBLE, 0, MPIWorld);
	}

	// INITIALIZATION: SAVING MESHES //

	Matrix* p_even     = new Matrix[(mycountX+2)*N_E];
	Matrix* m_even     = new Matrix[(mycountX+2)*N_E];
	Matrix* antip_even = new Matrix[(mycountX+2)*N_E];
	Matrix* antim_even = new Matrix[(mycountX+2)*N_E];

	Matrix* p_odd      = new Matrix[(mycountX+2)*N_E];
	Matrix* m_odd      = new Matrix[(mycountX+2)*N_E];
	Matrix* antip_odd  = new Matrix[(mycountX+2)*N_E];
	Matrix* antim_odd  = new Matrix[(mycountX+2)*N_E];

	// These arrays will become 'actual' ones in RK main cycle

	Matrix* p; Matrix* antip;
	Matrix* m; Matrix* antim;

	// and these will become 'old' ones

	Matrix* old_p; Matrix* old_antip;
	Matrix* old_m; Matrix* old_antim;

	for (int x = 0; x < mycountX + 2; ++x)
	{
		for (int e = 0; e < N_E; ++e)
		{
			p_odd[x*N_E + e].unpack(scatt_buffer + x*4*4*N_E + e*4*4 + 0*4);
			m_odd[x*N_E + e].unpack(scatt_buffer + x*4*4*N_E + e*4*4 + 1*4);

			antip_odd[x*N_E + e].unpack(scatt_buffer + x*4*4*N_E + e*4*4 + 2*4);
			antim_odd[x*N_E + e].unpack(scatt_buffer + x*4*4*N_E + e*4*4 + 3*4);
		}
	}

	// Spectra buffers; used in integrals

	double* Spec_e     = new double[N_E];
	double* Spec_x     = new double[N_E];

	double* Spec_antie = new double[N_E];
	double* Spec_antix = new double[N_E];

	double* Spec_nu     = new double[N_E];
	double* Spec_antinu = new double[N_E];

	// Buffer which contains the constant terms in Hamiltonians,
	// depends only on energy

	Matrix* H_0_mesh = new Matrix[N_E];

	// Saving spectra in a buffer, to reduce re-calculations for these values
	// and also initialize the constant terms of Hamiltonian

	for (int s = 0; s < N_E; ++s)
	{
		// current point of energy in spectra
		double E = Emin + s*dE;

		if (N_E != 1)
		{
			Spec_e[s]     = Sp(E, 0);
			Spec_x[s]     = Sp(E, 1);

			Spec_antie[s] = Sp(E, 2);
			Spec_antix[s] = Sp(E, 3);
		}
		else
		{
			Spec_e[s]     = alpha / (1.0 + alpha) * ee0     / dE;
			Spec_x[s]     = alpha / (1.0 + alpha) * xx0     / dE;

			Spec_antie[s] = 1.0   / (1.0 + alpha) * antiee0 / dE;
			Spec_antix[s] = 1.0   / (1.0 + alpha) * antixx0 / dE;
		}

		Spec_nu[s]     = Spec_e[s]     + Spec_x[s];
		Spec_antinu[s] = Spec_antie[s] + Spec_antix[s];

		// the part of Hamiltonian depending on E
		H_0_mesh[s] = -eta*dm2*c2theta/(4.0*E) * tau3 + eta*dm2*s2theta/(4.0*E) * tau1;
	}

	// CALCULATIONAL PART: RK-CYCLE //

	Matrix* K1 = new Matrix[(mycountX+2)*N_E * 4]; // {all p, m, antip, antim are put here, in one array}
	Matrix* K2 = new Matrix[(mycountX+2)*N_E * 4];
	Matrix* K3 = new Matrix[(mycountX+2)*N_E * 4];
	Matrix* K4 = new Matrix[(mycountX+2)*N_E * 4];

	// Temporary buffers used for interchanges

	Complex* leftsend  = new Complex[4 * N_E*4*4]; // {each nth N_E*4*4-block corresponds K_n coefficients}
	Complex* rightsend = new Complex[4 * N_E*4*4];

	Complex* leftrecv  = new Complex[4 * N_E*4*4];
	Complex* rightrecv = new Complex[4 * N_E*4*4];

	double time = 0;
	double time_here = 0;

	// Initialize variable of initial time stamp

	if (rank == 0) {
		time = MPI_Wtime();
	}

	// Initialize displacements and counts for MPI_Gather section

	get_line_counts(size, counts);
	get_line_displacements(size, displacements);

	// Initialize grid-files on master-node
	// and bin-files, if it's the beginning of caclulations, Z_init = 0

	if (rank == 0)
	{
		make_gridx();
		make_gridE();

		if (Z_init == 0)
		{
			make_bin_files();
			make_gridz();
		}
	}

	double Z;
	int ITER;

	for (int z = 0; z < N_Z; ++z)
	{
		for (int zi = 0; zi < STEP_Z; ++zi)
		{
			// Current index of iteration and Z-coordinate

			ITER = z * STEP_Z + zi;
			Z    = Z_init + (ITER + 1)*dZ;

			// Change pointers from old ones to buffers for writing
			// That's made for solving race-conditions in tick-RK4 section
			// while writing new data with using data of meshes

			if (ITER % 2 == 0)
			{
				p = p_even; antip = antip_even;
				m = m_even; antim = antim_even;

				old_p = p_odd; old_antip = antip_odd;
				old_m = m_odd; old_antim = antim_odd;
			}
			else
			{
				p = p_odd; antip = antip_odd;
				m = m_odd; antim = antim_odd;

				old_p = p_even; old_antip = antip_even;
				old_m = m_even; old_antim = antim_even;
			}

			// 1ST COEFFICIENT OF RK //

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

					K1[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s]);     // p
					K1[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s]);     // m
					K1[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s]); // antip
					K1[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s]); // antim

					// If it's the boundary of chunck, data will be put in the corresponding buffer
					// for sending to neighbour rank

					if (x == 1)
					{
						K1[4*x*N_E + 4*s + 0].pack(leftsend + 16*s + 0*4);
						K1[4*x*N_E + 4*s + 1].pack(leftsend + 16*s + 1*4);
						K1[4*x*N_E + 4*s + 2].pack(leftsend + 16*s + 2*4);
						K1[4*x*N_E + 4*s + 3].pack(leftsend + 16*s + 3*4);
					}
					else if (x == mycountX)
					{
						K1[4*x*N_E + 4*s + 0].pack(rightsend + 16*s + 0*4);
						K1[4*x*N_E + 4*s + 1].pack(rightsend + 16*s + 1*4);
						K1[4*x*N_E + 4*s + 2].pack(rightsend + 16*s + 2*4);
						K1[4*x*N_E + 4*s + 3].pack(rightsend + 16*s + 3*4);
					}
				}
			}

			// Interchanging the first coefficients of RK-cycle

			MPI_Request left_req_1, right_req_1;
			MPI_Status  left_st_1,  right_st_1;

			MPI_Isend(leftsend,  1, MPISpectra, left_neighbour,  1, MPIWorld,  &left_req_1);
			MPI_Isend(rightsend, 1, MPISpectra, right_neighbour, 1, MPIWorld, &right_req_1);

			MPI_Recv(leftrecv,   1, MPISpectra, left_neighbour,  1, MPIWorld,  &left_st_1);
			MPI_Recv(rightrecv,  1, MPISpectra, right_neighbour, 1, MPIWorld, &right_st_1);

			// Unpacking and saving on meshes 'K1'

			for (int s = 0; s < N_E; ++s)
			{
				int x = 0;

				K1[4*x*N_E + 4*s + 0].unpack(leftrecv + 16*s + 0*4);
				K1[4*x*N_E + 4*s + 1].unpack(leftrecv + 16*s + 1*4);
				K1[4*x*N_E + 4*s + 2].unpack(leftrecv + 16*s + 2*4);
				K1[4*x*N_E + 4*s + 3].unpack(leftrecv + 16*s + 3*4);

				x = mycountX + 1;

				K1[4*x*N_E + 4*s + 0].unpack(rightrecv + 16*s + 0*4);
				K1[4*x*N_E + 4*s + 1].unpack(rightrecv + 16*s + 1*4);
				K1[4*x*N_E + 4*s + 2].unpack(rightrecv + 16*s + 2*4);
				K1[4*x*N_E + 4*s + 3].unpack(rightrecv + 16*s + 3*4);	
			}

			// 2ND COEFFICIENT OF RK //

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

					K2[4*x*N_E + 4*s + 0] = -1.0*tanomega*Dp     - I/(ph::hbar*ph::c*comega)*comm(H_p,     old_p[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 0]*dZ);     // p
					K2[4*x*N_E + 4*s + 1] =      tanomega*Dm     - I/(ph::hbar*ph::c*comega)*comm(H_m,     old_m[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 1]*dZ);     // m
					K2[4*x*N_E + 4*s + 2] = -1.0*tanomega*Dantip - I/(ph::hbar*ph::c*comega)*comm(H_antip, old_antip[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 2]*dZ); // antip
					K2[4*x*N_E + 4*s + 3] =      tanomega*Dantim - I/(ph::hbar*ph::c*comega)*comm(H_antim, old_antim[x*N_E + s] + 0.5*K1[4*x*N_E + 4*s + 3]*dZ); // antim

					// If it's the boundary of chunck, data will be put in the corresponding buffer
					// for sending to neighbour rank

					if (x == 1)
					{
						K2[4*x*N_E + 4*s + 0].pack(leftsend + 16*N_E + 16*s + 0*4);
						K2[4*x*N_E + 4*s + 1].pack(leftsend + 16*N_E + 16*s + 1*4);
						K2[4*x*N_E + 4*s + 2].pack(leftsend + 16*N_E + 16*s + 2*4);
						K2[4*x*N_E + 4*s + 3].pack(leftsend + 16*N_E + 16*s + 3*4);
					}
					else if (x == mycountX)
					{
						K2[4*x*N_E + 4*s + 0].pack(rightsend + 16*N_E + 16*s + 0*4);
						K2[4*x*N_E + 4*s + 1].pack(rightsend + 16*N_E + 16*s + 1*4);
						K2[4*x*N_E + 4*s + 2].pack(rightsend + 16*N_E + 16*s + 2*4);
						K2[4*x*N_E + 4*s + 3].pack(rightsend + 16*N_E + 16*s + 3*4);
					}
				}
			}

			// Interchanging the second coefficients of RK-cycle

			MPI_Request left_req_2, right_req_2;
			MPI_Status  left_st_2,  right_st_2;

			MPI_Isend(leftsend  + 16*N_E, 1, MPISpectra, left_neighbour,  2, MPIWorld,  &left_req_2);
			MPI_Isend(rightsend + 16*N_E, 1, MPISpectra, right_neighbour, 2, MPIWorld, &right_req_2);

			MPI_Recv(leftrecv   + 16*N_E, 1, MPISpectra, left_neighbour,  2, MPIWorld,  &left_st_2);
			MPI_Recv(rightrecv  + 16*N_E, 1, MPISpectra, right_neighbour, 2, MPIWorld, &right_st_2);

			// Unpacking and saving on meshes 'K2'

			for (int s = 0; s < N_E; ++s)
			{
				int x = 0;

				K2[4*x*N_E + 4*s + 0].unpack(leftrecv + 16*N_E + 16*s + 0*4);
				K2[4*x*N_E + 4*s + 1].unpack(leftrecv + 16*N_E + 16*s + 1*4);
				K2[4*x*N_E + 4*s + 2].unpack(leftrecv + 16*N_E + 16*s + 2*4);
				K2[4*x*N_E + 4*s + 3].unpack(leftrecv + 16*N_E + 16*s + 3*4);

				x = mycountX + 1;

				K2[4*x*N_E + 4*s + 0].unpack(rightrecv + 16*N_E + 16*s + 0*4);
				K2[4*x*N_E + 4*s + 1].unpack(rightrecv + 16*N_E + 16*s + 1*4);
				K2[4*x*N_E + 4*s + 2].unpack(rightrecv + 16*N_E + 16*s + 2*4);
				K2[4*x*N_E + 4*s + 3].unpack(rightrecv + 16*N_E + 16*s + 3*4);	
			}

			// 3RD COEFFICIENT OF RK //

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

					// If it's the boundary of chunck, data will be put in the corresponding buffer
					// for sending to neighbour rank

					if (x == 1)
					{
						K3[4*x*N_E + 4*s + 0].pack(leftsend + 32*N_E + 16*s + 0*4);
						K3[4*x*N_E + 4*s + 1].pack(leftsend + 32*N_E + 16*s + 1*4);
						K3[4*x*N_E + 4*s + 2].pack(leftsend + 32*N_E + 16*s + 2*4);
						K3[4*x*N_E + 4*s + 3].pack(leftsend + 32*N_E + 16*s + 3*4);
					}
					else if (x == mycountX)
					{
						K3[4*x*N_E + 4*s + 0].pack(rightsend + 32*N_E + 16*s + 0*4);
						K3[4*x*N_E + 4*s + 1].pack(rightsend + 32*N_E + 16*s + 1*4);
						K3[4*x*N_E + 4*s + 2].pack(rightsend + 32*N_E + 16*s + 2*4);
						K3[4*x*N_E + 4*s + 3].pack(rightsend + 32*N_E + 16*s + 3*4);
					}
				}
			}

			// Interchanging the third coefficients of RK-cycle

			MPI_Request left_req_3, right_req_3;
			MPI_Status  left_st_3,  right_st_3;

			MPI_Isend(leftsend  + 32*N_E, 1, MPISpectra, left_neighbour,  3, MPIWorld,  &left_req_3);
			MPI_Isend(rightsend + 32*N_E, 1, MPISpectra, right_neighbour, 3, MPIWorld, &right_req_3);

			MPI_Recv(leftrecv   + 32*N_E, 1, MPISpectra, left_neighbour,  3, MPIWorld,  &left_st_3);
			MPI_Recv(rightrecv  + 32*N_E, 1, MPISpectra, right_neighbour, 3, MPIWorld, &right_st_3);

			// Unpacking and saving on meshes 'K3'

			for (int s = 0; s < N_E; ++s)
			{
				int x = 0;

				K3[4*x*N_E + 4*s + 0].unpack(leftrecv + 32*N_E + 16*s + 0*4);
				K3[4*x*N_E + 4*s + 1].unpack(leftrecv + 32*N_E + 16*s + 1*4);
				K3[4*x*N_E + 4*s + 2].unpack(leftrecv + 32*N_E + 16*s + 2*4);
				K3[4*x*N_E + 4*s + 3].unpack(leftrecv + 32*N_E + 16*s + 3*4);

				x = mycountX + 1;

				K3[4*x*N_E + 4*s + 0].unpack(rightrecv + 32*N_E + 16*s + 0*4);
				K3[4*x*N_E + 4*s + 1].unpack(rightrecv + 32*N_E + 16*s + 1*4);
				K3[4*x*N_E + 4*s + 2].unpack(rightrecv + 32*N_E + 16*s + 2*4);
				K3[4*x*N_E + 4*s + 3].unpack(rightrecv + 32*N_E + 16*s + 3*4);	
			}

			// 4TH COEFFICIENT OF RK //

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

					// If it's the boundary of chunck, data will be put in the corresponding buffer
					// for sending to neighbour rank

					if (x == 1)
					{
						K4[4*x*N_E + 4*s + 0].pack(leftsend + 48*N_E + 16*s + 0*4);
						K4[4*x*N_E + 4*s + 1].pack(leftsend + 48*N_E + 16*s + 1*4);
						K4[4*x*N_E + 4*s + 2].pack(leftsend + 48*N_E + 16*s + 2*4);
						K4[4*x*N_E + 4*s + 3].pack(leftsend + 48*N_E + 16*s + 3*4);
					}
					else if (x == mycountX)
					{
						K4[4*x*N_E + 4*s + 0].pack(rightsend + 48*N_E + 16*s + 0*4);
						K4[4*x*N_E + 4*s + 1].pack(rightsend + 48*N_E + 16*s + 1*4);
						K4[4*x*N_E + 4*s + 2].pack(rightsend + 48*N_E + 16*s + 2*4);
						K4[4*x*N_E + 4*s + 3].pack(rightsend + 48*N_E + 16*s + 3*4);
					}
				}
			}

			// Interchanging the fourth coefficients of RK-cycle

			MPI_Request left_req_4, right_req_4;
			MPI_Status  left_st_4,  right_st_4;

			MPI_Isend(leftsend  + 48*N_E, 1, MPISpectra, left_neighbour,  4, MPIWorld,  &left_req_4);
			MPI_Isend(rightsend + 48*N_E, 1, MPISpectra, right_neighbour, 4, MPIWorld, &right_req_4);

			MPI_Recv(leftrecv   + 48*N_E, 1, MPISpectra, left_neighbour,  4, MPIWorld,  &left_st_4);
			MPI_Recv(rightrecv  + 48*N_E, 1, MPISpectra, right_neighbour, 4, MPIWorld, &right_st_4);

			// Unpacking and saving on meshes 'K4'

			for (int s = 0; s < N_E; ++s)
			{
				int x = 0;

				K4[4*x*N_E + 4*s + 0].unpack(leftrecv + 48*N_E + 16*s + 0*4);
				K4[4*x*N_E + 4*s + 1].unpack(leftrecv + 48*N_E + 16*s + 1*4);
				K4[4*x*N_E + 4*s + 2].unpack(leftrecv + 48*N_E + 16*s + 2*4);
				K4[4*x*N_E + 4*s + 3].unpack(leftrecv + 48*N_E + 16*s + 3*4);

				x = mycountX + 1;

				K4[4*x*N_E + 4*s + 0].unpack(rightrecv + 48*N_E + 16*s + 0*4);
				K4[4*x*N_E + 4*s + 1].unpack(rightrecv + 48*N_E + 16*s + 1*4);
				K4[4*x*N_E + 4*s + 2].unpack(rightrecv + 48*N_E + 16*s + 2*4);
				K4[4*x*N_E + 4*s + 3].unpack(rightrecv + 48*N_E + 16*s + 3*4);	
			}

			// SAVING NEW LINE ON LOCAL MESH //

			for (int x = 0; x < mycountX+2; ++x)
			{
				for (int s = 0; s < N_E; ++s)
				{
					p[x*N_E + s] = old_p[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 0] + 2.0*K2[4*x*N_E + 4*s + 0] + 2.0*K3[4*x*N_E + 4*s + 0] + K4[4*x*N_E + 4*s + 0])*dZ;
					m[x*N_E + s] = old_m[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 1] + 2.0*K2[4*x*N_E + 4*s + 1] + 2.0*K3[4*x*N_E + 4*s + 1] + K4[4*x*N_E + 4*s + 1])*dZ;

					antip[x*N_E + s] = old_antip[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 2] + 2.0*K2[4*x*N_E + 4*s + 2] + 2.0*K3[4*x*N_E + 4*s + 2] + K4[4*x*N_E + 4*s + 2])*dZ;
					antim[x*N_E + s] = old_antim[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 3] + 2.0*K2[4*x*N_E + 4*s + 3] + 2.0*K3[4*x*N_E + 4*s + 3] + K4[4*x*N_E + 4*s + 3])*dZ;

					if (NORM)
					{
						// Recover hermitance and trace=1

						p[x*N_E + s].normalize();
						m[x*N_E + s].normalize();
						antip[x*N_E + s].normalize();
						antim[x*N_E + s].normalize();
					}
				}
			}
		}

		// GATHER NEW LINE ON MASTER NODE AND SAVE IT //

		// Copying data from arrays of Matrix to arrays of Complex

		for (int x = 1; x <= mycountX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				p[x*N_E + e].pack(gath_buffer     + (x-1)*N_E*16 + e*16 + 0*4);
				m[x*N_E + e].pack(gath_buffer     + (x-1)*N_E*16 + e*16 + 1*4);
				antip[x*N_E + e].pack(gath_buffer + (x-1)*N_E*16 + e*16 + 2*4);
				antim[x*N_E + e].pack(gath_buffer + (x-1)*N_E*16 + e*16 + 3*4);
			}
		}

		// Gathering and saving

		if (rank == 0)
		{
			MPI_Gatherv(MPI_IN_PLACE, mycountX * N_E * 16, MPIComplex, gath_buffer, counts, displacements, MPIComplex, 0, MPIWorld);

			Bin_output_line(size, counts, displacements, gath_buffer);
			Bin_output_rec(size, counts, displacements, gath_buffer);
			Text_output_Z_final(Z);

			add_to_gridz(Z);
		}
		else
		{
			MPI_Gatherv(gath_buffer, mycountX * N_E * 16, MPIComplex, NULL, counts, displacements, MPIComplex, 0, MPIWorld);			
		}

		// LOG THIS IN STDOUT //

		if (rank == 0)
		{
			time_here = MPI_Wtime();
			std::cout << std::fixed << std::setprecision(2) << "z=" << Z/ph::km << " of "
			          << (Z_init + Z_displacement)/ph::km << "; t=" << time_here - time << " sec" << std::endl;			
		}
	}
}