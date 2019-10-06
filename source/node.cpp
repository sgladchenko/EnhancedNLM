#include "node.h"
#include "segmentation.h"
#include "io.h"
#include "matrix.h"
#include "pmatrix.h"
#include "InitialSpectra.h"
#include "inhomogeneities.h"
#include "evolution.h"
#include "pmatrix.h"

#include <iostream>
#include <iomanip>

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

	// counts/displacements for gathering the adiabaticities

	int* ad_counts        = new int[size];
	int* ad_displacements = new int[size];

	get_ad_counts(size, ad_counts);
	get_ad_displacements(size, ad_displacements);

	// Scattering of the initial data

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

	// NEW: Grids of the Hamiltonians -- we're to estimate the adiabaticity of only one flow

	PMatrix* PH_p_even = new PMatrix[(mycountX + 1)*N_E]; // we actually need just the left neighbour's element
	PMatrix* PH_p_odd  = new PMatrix[(mycountX + 1)*N_E]; // due to we're employing the 1st order exressions

	PMatrix* PH_p_old;
	PMatrix* PH_p;

	// Buffers of gathering the adiabaticity factor

	Complex* adiabaticity_buffer;

	if (rank == 0)
	{
		adiabaticity_buffer = new Complex[N_X * N_E];		
	}
	else
	{
		adiabaticity_buffer = new Complex[mycountX * N_E];
	}

	// Extra buffers for interchanges of the Hamiltonians: in order to evaluate the adiabaticity factor

	Complex* send_PH_p_buf = new Complex[N_E * 4];
	Complex* recv_PH_p_buf = new Complex[N_E * 4];

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

	// Time stamps

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
		for (int zi = 0; zi < STEP_Z; ++zi) // The calculational phase
		{
			// Current index of iteration and Z-coordinate

			ITER = z * STEP_Z + zi;
			Z    = Z_init + (ITER + 1)*dZ;

			// Change pointers from old ones to buffers for writing
			// That's made for solving race-conditions in tick-RK4 section
			// while writing new data with using data of meshes

			if (ITER % 2 == 0)
			{
				p         = p_even;
				m         = m_even;
				antip     = antip_even;
				antim     = antim_even;

				old_p     = p_odd; 
				old_m     = m_odd;
				old_antip = antip_odd;
				old_antim = antim_odd;

				PH_p      = PH_p_even;
				PH_p_old  = PH_p_odd;
			}
			else
			{
				p         = p_odd; 
				m         = m_odd;
				antip     = antip_odd;
				antim     = antim_odd;

				old_p     = p_even; 
				old_m     = m_even; 
				old_antip = antip_even;
				old_antim = antim_even;

				PH_p      = PH_p_odd;
				PH_p_old  = PH_p_even;
			}

			// 1ST COEFFICIENT OF RK //

			evaluate_K1(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        harmonics,
				        mycountX,
				        Z,
				        Xleft,
				        PH_p,
				        K1);

			interchange_Ks(K1,
						   leftsend,
					       rightsend,
					       leftrecv,
					       rightrecv,
					       left_neighbour,
					       right_neighbour,
					       1,
					       mycountX);

			if (AD && ((zi == STEP_Z - 1) || (zi == STEP_Z - 2))) // Either it's the last iteration of this phase or before the last one
			{
				interchange_PH_p(PH_p,
				                 send_PH_p_buf,
				   			     recv_PH_p_buf,
					  			 left_neighbour,
					  			 right_neighbour,
					  			 mycountX);

				if (zi == STEP_Z - 1) // At the end of this calculational phase -- calculating
				{
					evaluate_adiabaticity(PH_p_old, PH_p, adiabaticity_buffer, mycountX);
				}
				else if (zi == STEP_Z - 2) // Before the last one -- something has to be prepared
				{
					// Explicit normalisation; it won't be conducted at this iteration
					// because evaluate_adiabaticity won't work there.
					// But at the next iteration we need exactly normalised values
					// in the PH_p_old array

					for (int x = 0; x <= mycountX; ++x)
					{
						for (int s = 0; s < N_E; ++s)
						{
							PH_p[x*N_E + s].normalise();
						}
					}
				}
			}

			// 2ND COEFFICIENT OF RK //

			evaluate_K2(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        harmonics,
				        mycountX,
				        Z,
				        Xleft,
				        K1,
				        K2);

			interchange_Ks(K2,
						   leftsend,
					       rightsend,
					       leftrecv,
					       rightrecv,
					       left_neighbour,
					       right_neighbour,
					       2,
					       mycountX);

			// 3RD COEFFICIENT OF RK //

			evaluate_K3(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        harmonics,
				        mycountX,
				        Z,
				        Xleft,
				        K2,
				        K3);

			interchange_Ks(K3,
						   leftsend,
					       rightsend,
					       leftrecv,
					       rightrecv,
					       left_neighbour,
					       right_neighbour,
					       3,
					       mycountX);

			// 4TH COEFFICIENT OF RK //

			evaluate_K4(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        harmonics,
				        mycountX,
				        Z,
				        Xleft,
				        K3,
				        K4);

			interchange_Ks(K4,
						   leftsend,
					       rightsend,
					       leftrecv,
					       rightrecv,
					       left_neighbour,
					       right_neighbour,
					       4,
					       mycountX);

			// SAVING NEW LINE IN THE LOCAL MESH //

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

		// GATHER NEW LINE AT THE MASTER NODE AND SAVE //

		// Copying data from the arrays of Matrix to the arrays of Complex

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

		// Gathering and saving the main calculational data

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

		// Gathering the adiabaticity grid, if needed

		if (AD)
		{
			if (rank == 0)
			{
				MPI_Gatherv(MPI_IN_PLACE, mycountX*N_E, MPIComplex, adiabaticity_buffer, ad_counts, ad_displacements, MPIComplex, 0, MPIWorld);

				// Saving at the binary file
				Bin_output_ad(size, ad_counts, ad_displacements, adiabaticity_buffer);
			}
			else
			{
				MPI_Gatherv(adiabaticity_buffer, mycountX*N_E, MPIComplex, NULL, ad_counts, ad_displacements, MPIComplex, 0, MPIWorld);
			}
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