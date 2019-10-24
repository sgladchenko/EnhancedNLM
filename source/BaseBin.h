#ifndef __BASEBIN_H
#define __BASEBIN_H

#include "Constants.h"

class BaseBin
// The base class of processing of the binary files with the basic functionality
// (i.e. there might some other classes based on this one with extended methods,
// e.g. to save some other information as the adiabaticity grid etc.)
{
	public:
		//BaseBin();
		BaseBin(int given_CommSize);
		~BaseBin();

		// Different inputs
		void InputScatterBuffer(Complex* ScatterBuffer);
		void InputZ_init(double& Z_init);
		void InputHarmonics(double* Harmonics);

		// Different outputs
		void OutputLine(Complex* GatherBuffer);
		void OutputRec(Complex* GatherBuffer);
		void OutputZ_final(double Z_final);

		// Generating files of the grids
		void MakeZgrid();
		void MakeXgrid();
		void MakeEgrid();
		void MakeBins();

		// Add a new point to the Z-grid file
		void AddToZgrid(double Z_here);

		// Size of the communicator
		int CommSize;

	protected:

		// Some other handy stuff
		int* CountXs;
		int* GatherCounts;
		int* GatherDisplacements;
};

// Filenames of the binary input/output

#define P_BIN "./data/p.bin"
#define M_BIN "./data/m.bin"
#define ANTIP_BIN "./data/antip.bin"
#define ANTIM_BIN "./data/antim.bin"

#define REC_BIN "./data/rec.bin"

#define ZGRID "./data/zgrid.bin"
#define XGRID "./data/xgrid.bin"
#define EGRID "./data/egrid.bin"

#define Z_INIT "./data/z.txt"

#define HARMONICS_BIN "./data/harmonics.bin"

#endif