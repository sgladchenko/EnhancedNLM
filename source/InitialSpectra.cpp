#include "InitialSpectra.h"


// The Fermi distribution f(x) = x^2 / (e^x + 1) achieves the maximum at x = xMax and f(xMax) = fMax
const Real InitialSpectra::xMax = 2.21771510575709;
const Real InitialSpectra::fMax = 0.48283007878900;

// Integral of f(x)dx from 0 to +infinity
const Real InitialSpectra::FermiCurveIntegral_dx_Reals = 1.803279;


InitialSpectra::InitialSpectra(Real _totalLuminosity)
    : totalLuminosity(_totalLuminosity)
{
	SetDefaultTemperatures();
	Normalize();
}


void InitialSpectra::Normalize()
{
	// We should have sum(f = 1..2*Nf, integral(s(E,f)dE)) = totalLuminosity
	Real norm = 0.0;
	for(unsigned f = 0; f < E_at_max.size(); f++)
	{
		norm += FermiCurveIntegral_dE_Reals(E_at_max[f], value_at_max[f]);
	}
	
	Real factor = totalLuminosity / norm;
	for(unsigned f = 0; f < E_at_max.size(); f++)
	{
		value_at_max[f] *= factor;
	}
}


void InitialSpectra::SetDefaultTemperatures()
{
    // The data taken from Fig.1 in [A. de Gouvea, S. Shalgar, Effect of Transition Magnetic Moments on Collective Supernova Neutrino Oscillations, 
	// arXiv:1207.0516 ]
	
	using PhysicalConstants::MeV;
	
	// nu_e
	E_at_max.push_back(8.0 * MeV);       
	value_at_max.push_back(3.35);    // arb. units
	// nu_x
	E_at_max.push_back(12.5 * MeV);       
	value_at_max.push_back(1.9);    // arb. units
	// antinu_e
	E_at_max.push_back(10.0 * MeV);       
	value_at_max.push_back(1.7);    // arb. units
	// antinu_x
	E_at_max.push_back(12.5 * MeV);       
	value_at_max.push_back(1.9);    // arb. units
}


/* static */ 
Real InitialSpectra::FermiCurveIntegral_dx(Real x1, Real x2)
{
    if(x1 == x2) return 0;
    if(x2 < x1) return -FermiCurveIntegral_dx(x2, x1);
	
    const Real dx = 0.01;
	Real integral = 0;
	for(Real x = x1; x < x2; x += dx)
	{
	    integral += FermiCurve(x) * dx;
	}
	
	return integral;
}


// The Fermi (Planck) distribution of the particle number over the neutrino energy 
// E is the neutrino energy, E_at_max is the position of the peak, value_at_max is the peak height
/* static */
Real InitialSpectra::FermiSpectrum(Real E, Real E_at_max, Real value_at_max)
{	  
	return (value_at_max / fMax) * FermiCurve((E / E_at_max) * xMax);
}

// integral of FermiSpectrum over dE from 0 to +infinity
/* static */
Real InitialSpectra::FermiCurveIntegral_dE_Reals(Real E_at_max, Real value_at_max)
{
	// Int(s(E)dE) = < E = theta * x > = theta * Int(s(x * theta)dx) = (value_at_max / fMax) * theta * FermiCurveIntegral_Reals
	Real theta = E_at_max / xMax;
	return (value_at_max / fMax) * theta * FermiCurveIntegral_dx_Reals;
}

// what was called s_f(E) in our papers
Real InitialSpectra::Spectrum(Real E, int flavor) const
{
    if(flavor < 0 || flavor >= (int)E_at_max.size())
		return 0;
	return FermiSpectrum(E, E_at_max[flavor], value_at_max[flavor]);
}


#ifdef __INITIALSPECTRA_TEST

#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>

using namespace std;

// Convert real numbers to Mathematica notation (e.g. 6.02*^23 instead of 6.02E23)
string RealToMathematicaNumber(double x)
{
	char buffer[32];
	sprintf(buffer, "%lg", x);
	// Find letter 'e' or 'E' in the result
	char* p = buffer;
    while(*p && *p != 'e' && *p != 'E')
		p++;
	
	if(*p == 0)  // No exponent, we have come to the ending '\0'
		return string(buffer);
     
	// Replace the exponent symbol with '\0' => split the string into two!
	*p = 0;
    p++;
	return string(buffer) + "*^" + string(p);
}


int main(void)
{	
	using PhysicalConstants::MeV;
	
    ofstream os("data\\initialSpectra.test.m");
	InitialSpectra s(1e53);
	Real Emin = 0, Emax = 50 * MeV, dE = 0.5 * MeV;
	
	os  << "(* Initial spectra of a protoneutron star with total luminosity L = " << s.Luminosity() << " neutrinos / sec          *)" << endl
	    << "(* The elements of the array are in the form { E[MeV], s_nu_e(E) [erg^(-1)], s_nu_x(E), s_antinu_e(E), s_antinu_x(E)} *)" << endl << endl
		<< "{" << endl;
		
	os  << "    {" << endl 
	    << "    	" << RealToMathematicaNumber(50.0 * 1e5) << " (* radius, in cm *), " << endl 
		<< "    	{" << endl;
	for(Real E = Emin; E <= Emax; E += dE)
	{
		
		os << "     	   { " << fixed << setprecision(3) << E / MeV << ", ";
		for(unsigned f = 0; f < 4; f++)
		{
			os << scientific << RealToMathematicaNumber(s(E, f)); 
			if(f < 3) os << ", ";
		}
		os << " }";
		if(E + dE <= Emax) os << "," << endl;
	}
	os << endl 
	   << "        }" << endl
	   << "    }" << endl
	   << "}" << endl;
	
	return 0;
}

#endif