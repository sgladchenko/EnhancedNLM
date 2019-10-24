#include "BaseScheme.h"
#include "BaseBin.h"

#include "InitialSpectra.h"
#include "Segmentation.h"
#include "Constants.h"
#include "Matrix.h"

class StandardScheme: public BaseScheme
// The implementation of the numerical scheme with
// included adiabaticity
{
	public:
		StandardScheme(int given_MyRank, int given_CommSize);

		void Initialise();
		void Process();

		~StandardScheme();

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
		Matrix* p_even;
		Matrix* m_even;
		Matrix* antip_even;
		Matrix* antim_even;

		Matrix* p_odd;
		Matrix* m_odd;
		Matrix* antip_odd;
		Matrix* antim_odd;

		// These arrays will become 'actual' ones in RK main cycle
		Matrix* p; Matrix* antip;
		Matrix* m; Matrix* antim;

		// and these will become 'old' ones
		Matrix* old_p; Matrix* old_antip;
		Matrix* old_m; Matrix* old_antim;

		// Spectra buffers; used in integrals
		double* Spec_e;
		double* Spec_x;

		double* Spec_antie;
		double* Spec_antix;

		double* Spec_nu;
		double* Spec_antinu;

		// Buffer which contains the constant terms in Hamiltonians,
		// depends only on energy
		Matrix* H_0_mesh;

		// The very start time stamp of the calculations
		double Start_time;

		// Some temporary buffers used for interchanges
		Complex* LeftSend_buf;
		Complex* RightSend_buf;

		Complex* LeftRecv_buf;
		Complex* RightRecv_buf;

		// Matrices of the derivatives in the RK algorithm
		Matrix* K1;
		Matrix* K2;
		Matrix* K3;
		Matrix* K4;
};