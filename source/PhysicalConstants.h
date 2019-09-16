#ifndef __PHYSICALCONSTANTS_H
#define __PHYSICALCONSTANTS_H

#include <complex>
#include <math.h>
#include <iostream>

namespace PhysicalConstants
{
  typedef double Real;
  typedef std::complex<Real> Complex;

  const Complex imagUnit(0,1); 

  const Real Pi = 3.1415926535897932384626433832795;
  const Real Degree = Pi / 180;                      // in radians
  
  const Real cm   = 1.0;
  const Real sec  = 1.0;
  const Real gram = 1.0;
  const Real erg  = 1.0;
  const Real statCoulomb = 1.0;
  const Real Gauss = 1.0;
  
  const Real ElectronVolt = 1.602e-12 * erg;
  const Real eV = ElectronVolt;
  const Real MeV = 1e6 * eV;
  const Real FermiConstant = 1.166e-11 / (MeV*MeV); 
  const Real ProtonMass = 938.272 * MeV;
  const Real NeutronMass = 939.565 * MeV;
  const Real ElectronMass = 0.511 * MeV;
  const Real PlanckConstantReduced = 1.0546e-27 * erg * sec;
  const Real SpeedOfLight = 2.9979e+10 * cm / sec;
  const Real AvogadroConstant = 6.0221e23; // mol^(-1)
  const Real BohrMagneton = 9.27402e-21 * (statCoulomb * cm);
  const Real SchwingerField = 4.414e13 * Gauss;
  const Real WeakMixingAngle = asin(sqrt(0.24));
  
  // Aliases
  const Real GF = FermiConstant;
  const Real hbar = PlanckConstantReduced;
  const Real c = SpeedOfLight;
  const Real NA = AvogadroConstant;
  const Real muB = BohrMagneton;
  const Real H0 = SchwingerField;
  const Real km = 1e5 * cm;
  const Real m_e = ElectronMass;
  const Real theta_W = WeakMixingAngle;

  // Neutrino masses and mixing
  const Real deltaMSquared = 2.5e-3 * eV * eV; // 7.65e-5 * (eV * eV); // Mass-squared difference, in erg^2
  //const Real theta0 = asin(sqrt(0.304)); // Vacuum mixing angle, dimensionless
  const Real theta0 = 9 * Degree;   // Vacuum mixing angle, dimensionless

  
  // m1-m2 matrix element of the (Majorana) neutrino magnetic moment (tentative value)
  const Real mu12 = 1e-15 * muB;

  
  // Some easy functions
  inline Real sqr(Real x) { return x * x; }
  // Round to nearest integer
  inline int Round(Real x)
  {
     int xFloor = (int)floor(x);
     return x - xFloor < 0.5 ? xFloor : xFloor + 1;
  }
  // Arccotangent giving the result in (0, pi)
  inline Real acot(Real x) { return Pi/2 - atan(x); } 
  // Cotangent
  inline Real cot(Real x) { return tan(Pi/2 - x); }
  // Signum function
  inline int sgn(Real x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
  // The theta function
  inline int HeavisideTheta(Real x) { return x >= 0 ? 1 : 0; }

  // Write a complex number to a stream using "a+bi" notation
  inline void WriteComplex(std::ostream& os, Complex z)
  {
      if(norm(z) == 0) os << "0";
	  else if(z.imag() == 0) os << z.real();
	  else if(z.real() == 0) os << z.imag() << "i";
	  else os << z.real() << (z.imag() > 0 ? '+' : '-') << std::fabs(z.imag()) << "i";
  }

  inline Real amplitude(Complex c) 
  {
	  Real a, b, e;
	  a = c.real();
	  b = c.imag();
	  e = sqrt(a*a + b*b);
	  return e; 
  }

}

#endif