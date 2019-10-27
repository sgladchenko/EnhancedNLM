#include "AdiabaticScheme.h"
#include "Log.h"

AdiabaticScheme::AdiabaticScheme(int given_MyRank, int given_CommSize)
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
	H_0 = new Unit[N_E];

	// Define the first set of aliases
	Rho_prev   = Rhos;
	Rho_here   = Rhos + 1;
	Rho_next_0 = Rhos + 2;
	Rho_next_1 = Rhos + 3;

	H_prev   = Hs;
	H_here   = Hs + 1;
	H_next_0 = Hs + 2;
	H_next_1 = Hs + 3;

	// Allocation of the interchange buffers
	ToLeftNeighbour  = new Complex[16*N_E];
	ToRightNeighbour = new Complex[16*N_E];

	FromLeftNeighbour  = new Complex[16*N_E];
	FromRightNeighbour = new Complex[16*N_E];
}

void AdiabaticScheme::Initialise()
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
			Rho_prev->Init(MyCountX, ScatterBuffer);
			Rho_here->Init(MyCountX, ScatterBuffer);

			// An initial approximation of the next point at the 1st iteration
			Rho_next_0->Init(MyCountX, ScatterBuffer);
			Rho_next_1->Init(MyCountX, ScatterBuffer);

			// Remark: Layer class considers that all the grids should have additional points on the left and on the right
			//         from the neighbour-nodes

			// This is just a dummy intialisation -- in order to allocate desired memory for the further calculations

			H_prev->Init(MyCountX, ScatterBuffer);
			H_here->Init(MyCountX, ScatterBuffer);

			// An initial approximation of the next point at the 1st iteration
			H_next_0->Init(MyCountX, ScatterBuffer);
			H_next_1->Init(MyCountX, ScatterBuffer);
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
		H_0[s] = -eta*dm2*c2theta/(4.0*E) * sigma_3 + eta*dm2*s2theta/(4.0*E) * sigma_1; // Units
	}

	// Evaluating the initial Hamiltonians
	EvaluateHamiltonian(Rho_prev, H_prev, H_0,
			     		Spec_nu,
				  		Spec_antinu,
				  		Harmonics,
				  		MyCountX,
				  		0,
				  		Xleft);

	EvaluateHamiltonian(Rho_here, H_here, H_0,
			     		Spec_nu,
				  		Spec_antinu,
				  		Harmonics,
				  		MyCountX,
				  		0,
				  		Xleft);

	EvaluateHamiltonian(Rho_next_0, H_next_0, H_0,
			     		Spec_nu,
				  		Spec_antinu,
				  		Harmonics,
				  		MyCountX,
				  		0,
				  		Xleft);

	EvaluateHamiltonian(Rho_next_1, H_next_1, H_0,
			     		Spec_nu,
				  		Spec_antinu,
				  		Harmonics,
				  		MyCountX,
				  		0,
				  		Xleft);
}

void AdiabaticScheme::Process()
{
	// Object to operate the binary data
	BaseBin Bin(CommSize);

	// Logger object of the calculational process
	Log Logger(MyRank, Start_time, Z_init, Z_displacement);

	// Auxiliary variables
	double Z; int Step;

	for (int z = 0; z < N_Z; ++z)
	{
		// 'Calculational phase', and after each of them there's a gathering phase
		for (int zi = 0; zi < STEP_Z; ++zi)
		{
			// Current index of the step and Z-coordinate
			Step = z * STEP_Z + zi;
			Z    = Z_init + (Step)*dZ;

			// Calculations at this step

			AdiabaticFDM(Rho_prev, Rho_here, Rho_next_0, Rho_next_1,
						 H_prev,   H_here,   H_next_0,   H_next_1,
						 H_0,
						 Spec_nu,
						 Spec_antinu,
						 Harmonics,
						 MyCountX,
						 Z,
						 Xleft);

			// Pack new values to the interchange buffers
			(*Rho_next_0)(1).Pack(ToLeftNeighbour);
			(*Rho_next_0)(MyCountX).Pack(ToRightNeighbour);

			// Some stuff to be used in the interchanging processes
			MPI_Request left_req, right_req;
			MPI_Status  left_st,  right_st;

			// Interchanging...

			MPI_Isend(ToLeftNeighbour,  16*N_E, MPIComplex, LeftNeighbour,  1, MPIWorld,  &left_req);
			MPI_Isend(ToRightNeighbour, 16*N_E, MPIComplex, RightNeighbour, 1, MPIWorld, &right_req);

			MPI_Recv(FromLeftNeighbour,  16*N_E, MPIComplex, LeftNeighbour,  1, MPIWorld,  &left_st);
			MPI_Recv(FromRightNeighbour, 16*N_E, MPIComplex, RightNeighbour, 1, MPIWorld, &right_st);

			// Unpack new values for the interchange buffers
			(*Rho_next_0)(0).Unpack(FromLeftNeighbour);
			(*Rho_next_0)(MyCountX+1).Unpack(FromRightNeighbour);

			// Evaluate new Hamiltonians
			EvaluateHamiltonian(Rho_next_0, H_next_0, H_0,
					     		Spec_nu,
						  		Spec_antinu,
						  		Harmonics,
						  		MyCountX,
						  		Z,
						  		Xleft);

			// Cyclic swapping variables!
			Layer* tmp;

			tmp        = Rho_prev;
			Rho_prev   = Rho_here;
			Rho_here   = Rho_next_0;
			Rho_next_0 = tmp;

			tmp      = H_prev;
			H_prev   = H_here;
			H_here   = H_next_0;
			H_next_0 = tmp;
		}

		// Copying data from the arrays of Matrix to the arrays of Complex
		Rho_here->Dump(GatherBuffer);

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

AdiabaticScheme::~AdiabaticScheme()
{
	delete [] GatherBuffer;
	delete [] ScatterBuffer;
	delete [] Harmonics;

	delete [] Spec_e;
	delete [] Spec_x;

	delete [] Spec_antie;
	delete [] Spec_antix;

	delete [] Spec_nu;
	delete [] Spec_antinu;

	delete [] H_0;
}