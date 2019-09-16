#ifndef __INITIALSPECTRA_H
#define __INITIALSPECTRA_H

#include <vector>
#include <math.h>
#pragma hdrstop
#include "PhysicalConstants.h"

using PhysicalConstants::Real;

class InitialSpectra
{
    std::vector<Real> E_at_max, value_at_max;    // Real[2*Nf]
    Real totalLuminosity;                        // in particles / second
public:
    // Construct with default temperatures and ratios of flavor luminosities
	InitialSpectra(Real _totalLuminosity = 1.0);
	
	Real Luminosity() const { return totalLuminosity; }
	void SetLuminosity(Real newLumin) 
	{
		totalLuminosity = newLumin;
		Normalize();
	}
	
    // what was called s_f(E) in our papers
    Real Spectrum(Real E, int flavor) const;
	Real operator() (Real E, int flavor) const { return Spectrum(E, flavor); }
protected:
	void Normalize();
	void SetDefaultTemperatures();
	
	static Real FermiCurve(Real x)
	{
		// For safety 
		if(x > 100 || x < 0) return 0;
		return x * x / (exp(x) + 1);
	}
	
	// The Fermi (Planck) distribution of the particle number over the neutrino energy 
	// E is the neutrino energy, E_at_max is the position of the peak, value_at_max is the peak height
	static Real FermiSpectrum(Real E, Real E_at_max, Real value_at_max);
	
	// integral of FermiCurve from x1 to x2
	static Real FermiCurveIntegral_dx(Real x1, Real x2);
	
	// integral of FermiSpectrum over dE from 0 to +infinity
	static Real FermiCurveIntegral_dE_Reals(Real E_at_max, Real value_at_max);

private:
	// The Fermi distribution f(x) = x^2 / (e^x + 1) achieves the maximum at x = xMax and f(xMax) = fMax
	static const Real xMax, fMax;
	// Integral over dx from 0 to +infinity
	static const Real FermiCurveIntegral_dx_Reals;
};


#endif