#include "StandardScheme.h"
#include "StandardRK.h"
#include "Log.h"

StandardScheme::StandardScheme(int given_MyRank, int given_CommSize)
{
	MyRank = given_MyRank;
	CommSize = given_CommSize;

	// Determine my neighbour-nodes
	if (MyRank == 0)
	{
		LeftNeighbour  = CommSize - 1;
		RightNeighbour = 1;
	}
	else if (MyRank == CommSize - 1)
	{
		LeftNeighbour  = CommSize - 2;
		RightNeighbour = 0;
	}
	else
	{
		LeftNeighbour  = MyRank - 1;
		RightNeighbour = MyRank + 1;
	}

	// InitialSpectra object to be used in the calculations
	InitialSpectra Sp;

	// My chunck of the X-grid
	Obtain_SegmentX(MyRank, CommSize, MyLeftX, MyRightX);
	MyCountX = MyRightX - MyLeftX + 1;

	Xleft  = Xmin + MyLeftX*dX;
	Xright = Xmin + MyRightX*dX;

	// Allocations of the scattering and gathering buffers
	if (MyRank == 0)
	{
		GatherBuffer  = new Complex[N_X * 16*N_E];
		ScatterBuffer = new Complex[(N_X + 2*CommSize) * 16*N_E];
	}
	else
	{
		GatherBuffer  = new Complex[MyCountX * 16*N_E];
		ScatterBuffer = new Complex[(MyCountX + 2) * 16*N_E];
	}

	// Allocation of the Harmonics array
	if (NUM_HARMONICS > 0)
	{
		Harmonics = new double[NUM_HARMONICS];
	}

	// Arrays for the scattering processes
	ScatterCounts        = new int[CommSize];
	ScatterDisplacements = new int[CommSize];

	Obtain_ScatterCounts(CommSize, ScatterCounts);
	Obtain_ScatterDisplacements(CommSize, ScatterDisplacements);

	// Arrays for the gathering processes
	GatherCounts        = new int[CommSize];
	GatherDisplacements = new int[CommSize];

	Obtain_GatherCounts(CommSize, GatherCounts);
	Obtain_GatherDisplacements(CommSize, GatherDisplacements);

	// Data to be used in the calculations
	p_even     = new Matrix[(MyCountX+2)*N_E];
	m_even     = new Matrix[(MyCountX+2)*N_E];
	antip_even = new Matrix[(MyCountX+2)*N_E];
	antim_even = new Matrix[(MyCountX+2)*N_E];

	p_odd      = new Matrix[(MyCountX+2)*N_E];
	m_odd      = new Matrix[(MyCountX+2)*N_E];
	antip_odd  = new Matrix[(MyCountX+2)*N_E];
	antim_odd  = new Matrix[(MyCountX+2)*N_E];

	// Start point of the calculations
	Start_time = MPI_Wtime();

	// Spectra buffers; used in integrals
	Spec_e      = new double[N_E];
	Spec_x      = new double[N_E];

	Spec_antie  = new double[N_E];
	Spec_antix  = new double[N_E];

	Spec_nu     = new double[N_E];
	Spec_antinu = new double[N_E];

	// Buffer which contains the constant terms in Hamiltonians,
	// depends only on energy
	H_0_mesh = new Matrix[N_E];

	// Some temporary buffers used for interchanges
	LeftSend_buf  = new Complex[4 * N_E*4*4]; // {each nth N_E*4*4-block corresponds K_n coefficients}
	RightSend_buf = new Complex[4 * N_E*4*4];

	LeftRecv_buf  = new Complex[4 * N_E*4*4];
	RightRecv_buf = new Complex[4 * N_E*4*4];

	// Matrices of the derivatives in the RK algorithm
	K1 = new Matrix[(MyCountX+2)*N_E * 4]; // {all p, m, antip, antim are put here, in one array}
	K2 = new Matrix[(MyCountX+2)*N_E * 4];
	K3 = new Matrix[(MyCountX+2)*N_E * 4];
	K4 = new Matrix[(MyCountX+2)*N_E * 4];
}

void StandardScheme::Initialise()
{
	// Object to operate the binary data
	BaseBin Bin(CommSize);

	// Logger object of the initialisation process
	Log Logger(MyRank, Start_time, Z_init, Z_displacement);

	// Obtaining data from the binaries
	Logger.Out("Obtaining data from the binaries...");
	if (MyRank == 0)
	{
		Bin.InputScatterBuffer(ScatterBuffer);
		Bin.InputZ_init(Z_init);
		Bin.InputHarmonics(Harmonics);
	}

	// Initialize grid-files on master-node
	// and bin-files, if it's the beginning of the caclulations, when Z_init = 0
	Logger.Out("Generating binaries of the grids...");
	if (MyRank == 0)
	{
		Bin.MakeXgrid();
		Bin.MakeEgrid();

		if (Z_init == 0)
		{
			Bin.MakeBins();
			Bin.MakeZgrid();
		}
	}

	// Scattering the initial data
	Logger.Out("Scattering and casting the initial data...");
	MPI_Scatterv(ScatterBuffer,
				 ScatterCounts,
				 ScatterDisplacements,
				 MPIComplex,
				 ScatterBuffer,
				 (MyCountX + 2) * 16*N_E,
				 MPIComplex,
				 0,
				 MPIWorld);

	MPI_Bcast(&Z_init, 1, MPI_DOUBLE, 0, MPIWorld);
	MPI_Bcast(Harmonics, NUM_HARMONICS, MPI_DOUBLE, 0, MPIWorld);

	// Unpacking the obtained data from scattering
	Logger.Out("Unpacking the intial data...");
	for (int x = 0; x < MyCountX + 2; ++x)
	{
		for (int e = 0; e < N_E; ++e)
		{
			p_odd[x*N_E + e].unpack(ScatterBuffer + x*4*4*N_E + e*4*4 + 0*4);
			m_odd[x*N_E + e].unpack(ScatterBuffer + x*4*4*N_E + e*4*4 + 1*4);

			antip_odd[x*N_E + e].unpack(ScatterBuffer + x*4*4*N_E + e*4*4 + 2*4);
			antim_odd[x*N_E + e].unpack(ScatterBuffer + x*4*4*N_E + e*4*4 + 3*4);
		}
	}

	Logger.Out("Some other intialisations...");

	// Saving spectra in a buffer, to reduce re-calculations for these values
	// and also initialise the constant terms of Hamiltonian
	for (int s = 0; s < N_E; ++s)
	{
		// The current point of energy in spectra
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

		// The desired part of the Hamiltonian depending on E
		H_0_mesh[s] = -eta*dm2*c2theta/(4.0*E) * tau3 + eta*dm2*s2theta/(4.0*E) * tau1;
	}
}

void StandardScheme::Process()
{
	// Object to operate the binary data
	BaseBin Bin(CommSize);

	// Logger object of the calculational process
	Log Logger(MyRank, Start_time, Z_init, Z_displacement);

	// Auxiliary variables
	double Z; int ITER;

	for (int z = 0; z < N_Z; ++z)
	{
		for (int zi = 0; zi < STEP_Z; ++zi) // The 'calculational phase'
		{
			// Current index of iteration and Z-coordinate
			ITER = z * STEP_Z + zi;
			Z    = Z_init + (ITER + 1)*dZ;

			// Change pointers from old ones to buffers for writing
			// That's made for solving race-conditions in tick-RK4 section
			// while writing new data with using data of the grids
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
			}

			// Calculations at this step

			Evaluate_K1(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        Harmonics,
				        MyCountX,
				        Z,
				        Xleft,
				        K1);

			Interchange_Ks(K1,
						   LeftSend_buf,
					       RightSend_buf,
					       LeftRecv_buf,
					       RightRecv_buf,
					       LeftNeighbour,
					       RightNeighbour,
					       1,
					       MyCountX);

			Evaluate_K2(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        Harmonics,
				        MyCountX,
				        Z,
				        Xleft,
				        K1,
				        K2);

			Interchange_Ks(K2,
						   LeftSend_buf,
					       RightSend_buf,
					       LeftRecv_buf,
					       RightRecv_buf,
					       LeftNeighbour,
					       RightNeighbour,
					       2,
					       MyCountX);

			Evaluate_K3(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        Harmonics,
				        MyCountX,
				        Z,
				        Xleft,
				        K2,
				        K3);

			Interchange_Ks(K3,
						   LeftSend_buf,
					       RightSend_buf,
					       LeftRecv_buf,
					       RightRecv_buf,
					       LeftNeighbour,
					       RightNeighbour,
					       3,
					       MyCountX);

			Evaluate_K4(old_p,
 				 	    old_m,
				        old_antip,
				        old_antim,
				        H_0_mesh,
				        Spec_nu,
				        Spec_antinu,
				        Harmonics,
				        MyCountX,
				        Z,
				        Xleft,
				        K3,
				        K4);

			Interchange_Ks(K4,
						   LeftSend_buf,
					       RightSend_buf,
					       LeftRecv_buf,
					       RightRecv_buf,
					       LeftNeighbour,
					       RightNeighbour,
					       4,
					       MyCountX);

			// Calculate new points & Save new values in the buffers

			for (int x = 0; x < MyCountX+2; ++x)
			{
				for (int s = 0; s < N_E; ++s)
				{
					p[x*N_E + s] = old_p[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 0] + 2.0*K2[4*x*N_E + 4*s + 0] + 2.0*K3[4*x*N_E + 4*s + 0] + K4[4*x*N_E + 4*s + 0])*dZ;
					m[x*N_E + s] = old_m[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 1] + 2.0*K2[4*x*N_E + 4*s + 1] + 2.0*K3[4*x*N_E + 4*s + 1] + K4[4*x*N_E + 4*s + 1])*dZ;

					antip[x*N_E + s] = old_antip[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 2] + 2.0*K2[4*x*N_E + 4*s + 2] + 2.0*K3[4*x*N_E + 4*s + 2] + K4[4*x*N_E + 4*s + 2])*dZ;
					antim[x*N_E + s] = old_antim[x*N_E + s] + 1.0 / 6.0 * (K1[4*x*N_E + 4*s + 3] + 2.0*K2[4*x*N_E + 4*s + 3] + 2.0*K3[4*x*N_E + 4*s + 3] + K4[4*x*N_E + 4*s + 3])*dZ;

					if (NORM)
					{
						// Recover hermitance and trace=1, if needed
						p[x*N_E + s].normalise();
						m[x*N_E + s].normalise();
						antip[x*N_E + s].normalise();
						antim[x*N_E + s].normalise();
					}
				}
			}
		}

		// Gather new line of values at the root-node and save in the binaries

		// Copying data from the arrays of Matrix to the arrays of Complex
		for (int x = 1; x <= MyCountX; ++x)
		{
			for (int e = 0; e < N_E; ++e)
			{
				p[x*N_E + e].pack(GatherBuffer     + (x-1)*N_E*16 + e*16 + 0*4);
				m[x*N_E + e].pack(GatherBuffer     + (x-1)*N_E*16 + e*16 + 1*4);
				antip[x*N_E + e].pack(GatherBuffer + (x-1)*N_E*16 + e*16 + 2*4);
				antim[x*N_E + e].pack(GatherBuffer + (x-1)*N_E*16 + e*16 + 3*4);
			}
		}

		// Gathering and saving main calculational data
		if (MyRank == 0)
		{
			MPI_Gatherv(MPI_IN_PLACE, MyCountX * N_E * 16, MPIComplex, GatherBuffer, GatherCounts, GatherDisplacements, MPIComplex, 0, MPIWorld);

			// Update the binaries of the functions' grids
			Bin.OutputLine(GatherBuffer);
			Bin.OutputRec(GatherBuffer);
			Bin.OutputZ_final(Z);

			// Add a new point to the array of Z-points -- to indicate that we've got this point
			Bin.AddToZgrid(Z);
		}
		else
		{
			MPI_Gatherv(GatherBuffer, MyCountX * N_E * 16, MPIComplex, NULL, GatherCounts, GatherDisplacements, MPIComplex, 0, MPIWorld);			
		}

		// Log that we've got this point
		Logger.OutIter(Z);
	}
}

StandardScheme::~StandardScheme()
{
	delete [] GatherBuffer;
	delete [] ScatterBuffer;
	delete [] Harmonics;

	delete [] p_even;
	delete [] m_even;
	delete [] antip_even;
	delete [] antim_even;

	delete [] p_odd;
	delete [] m_odd;
	delete [] antip_odd;
	delete [] antim_odd;

	delete [] Spec_e;
	delete [] Spec_x;

	delete [] Spec_antie;
	delete [] Spec_antix;

	delete [] Spec_nu;
	delete [] Spec_antinu;

	delete [] H_0_mesh;

	delete [] LeftSend_buf;
	delete [] RightSend_buf;

	delete [] LeftRecv_buf;
	delete [] RightRecv_buf;
}