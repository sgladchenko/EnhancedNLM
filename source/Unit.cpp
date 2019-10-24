#include "Containers.h"

// Pauli matrices
Unit Id(1, 0, 0, 1);
Unit sigma_1(0, 1, 1, 0);
Unit sigma_2(0, -Complex(0, 1.0), Complex(0, 1.0), 0);
Unit sigma_3(1, 0, 0, -1);

Unit PauliMatrices[4] = {Id, sigma_1, sigma_2, sigma_3};

Unit::Unit(Complex* buffer)
{
	// Each of our matrices consists of four elememts (2x2)
	for (int i = 0; i < 4; ++i)
	{
		ComplexBuffer[i] = buffer[i];
	}
}

Unit::Unit(Complex e0, Complex e1, Complex e2, Complex e3)
{
	// Must be suitable for some tasks, e.g. when we're posing the Pauli matrices
	ComplexBuffer[0] = e0;
	ComplexBuffer[1] = e1;
	ComplexBuffer[2] = e2;
	ComplexBuffer[3] = e3;
}

Unit::Unit()
{
	for (int i = 0; i < 4; ++i)
	{
		ComplexBuffer[i] = 0;
	}	
}

Unit::Unit(const Unit& u)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] = u.ComplexBuffer[k];
	}
}

void Unit::Init(Complex e0, Complex e1, Complex e2, Complex e3)
{
	// Must be suitable for some tasks, e.g. after allocations
	ComplexBuffer[0] = e0;
	ComplexBuffer[1] = e1;
	ComplexBuffer[2] = e2;
	ComplexBuffer[3] = e3;	
}

void Unit::Init(Complex* buffer)
{
	for (int i = 0; i < 4; ++i)
	{
		ComplexBuffer[i] = buffer[i];
	}
}

Unit Unit::operator+(const Unit& u) const
{
	return Unit(ComplexBuffer[0] + u.ComplexBuffer[0],
				ComplexBuffer[1] + u.ComplexBuffer[1],
				ComplexBuffer[2] + u.ComplexBuffer[2],
				ComplexBuffer[3] + u.ComplexBuffer[3]);
}

Unit Unit::operator-(const Unit& u) const
{
	return Unit(ComplexBuffer[0] - u.ComplexBuffer[0],
				ComplexBuffer[1] - u.ComplexBuffer[1],
				ComplexBuffer[2] - u.ComplexBuffer[2],
				ComplexBuffer[3] - u.ComplexBuffer[3]);
}

Unit Unit::operator*(const Unit& u) const
{
	return Unit(ComplexBuffer[0]*u.ComplexBuffer[0] + ComplexBuffer[1]*u.ComplexBuffer[2],
				ComplexBuffer[0]*u.ComplexBuffer[1] + ComplexBuffer[1]*u.ComplexBuffer[3],
				ComplexBuffer[2]*u.ComplexBuffer[0] + ComplexBuffer[3]*u.ComplexBuffer[2],
				ComplexBuffer[2]*u.ComplexBuffer[1] + ComplexBuffer[3]*u.ComplexBuffer[3]);
}

Unit Unit::operator/(const Complex& z) const
{
	return Unit(ComplexBuffer[0] / z,
				ComplexBuffer[1] / z,
				ComplexBuffer[2] / z,
				ComplexBuffer[3] / z);
}

// Right multiplication on a Complex number
Unit Unit::operator*(const Complex z) const
{
	return Unit(ComplexBuffer[0] * z,
				ComplexBuffer[1] * z,
				ComplexBuffer[2] * z,
				ComplexBuffer[3] * z);
}

// Left multiplication on a Complex number
Unit operator*(const Complex& z, const Unit& u)
{
	return Unit(u.ComplexBuffer[0] * z,
				u.ComplexBuffer[1] * z,
				u.ComplexBuffer[2] * z,
				u.ComplexBuffer[3] * z);
}

void Unit::operator+=(const Unit& u)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] += u.ComplexBuffer[k];
	}
}

void Unit::operator-=(const Unit& u)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] -= u.ComplexBuffer[k];
	}
}

void Unit::operator*=(const Complex& z)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] *= z;
	}
}

void Unit::operator/=(const Complex& z)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] /= z;
	}	
}

void Unit::operator=(const Unit& u)
{
	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] = u.ComplexBuffer[k];
	}
}

Complex Unit::Tr() const
{
	return ComplexBuffer[0] + ComplexBuffer[3];
}

void Unit::normalise()
{
	// Recover Hermitance
	Complex o0, o1, o2, o3; // old ones

	o0 = ComplexBuffer[0];
	o1 = ComplexBuffer[1];
	o2 = ComplexBuffer[2];
	o3 = ComplexBuffer[3];

	ComplexBuffer[0] = 0.5*(o0 + std::conj(o0));
	ComplexBuffer[1] = 0.5*(o1 + std::conj(o2));
	ComplexBuffer[2] = 0.5*(o2 + std::conj(o1));
	ComplexBuffer[3] = 0.5*(o3 + std::conj(o3));

	// Recover trace = 1
	Complex tr = ComplexBuffer[0] + ComplexBuffer[3];

	for (int k = 0; k < 4; ++k)
	{
		ComplexBuffer[k] /= tr;
	}
}

Complex& Unit::operator()(int i, int j)
{
	return ComplexBuffer[i*2 + j];
}

Complex& Unit::operator[](int k)
{
	return ComplexBuffer[k];
}

Complex Unit::norm()
{
	Complex Norm = 0;

	for (int k = 0; k < 4; ++k)
	{
		Norm += std::norm(ComplexBuffer[k]);
	}

	return Norm;
}

// Commutator of two matrices and their normalised skew-product

Unit commutator(const Unit& u1, const Unit& u2)
{
	return u1*u2 - u2*u1;
}

Unit skew(const Unit& u1, const Unit& u2)
{
	// Define their components in the Pauli basis

	Complex PauliBuffer1[4];
	Complex PauliBuffer2[4];

	for (int k = 0; k < 4; ++k)
	{
		Unit tmp1 = u1*PauliMatrices[k];
		Unit tmp2 = u2*PauliMatrices[k];

		PauliBuffer1[k] = tmp1.Tr();
		PauliBuffer2[k] = tmp2.Tr();
	}

	// Calculate their 3D-norms

	Complex Norm3D_1 = 0;
	Complex Norm3D_2 = 0;

	for (int k = 1; k < 4; ++k)
	{
		Norm3D_1 += std::norm(PauliBuffer1[k]);
		Norm3D_2 += std::norm(PauliBuffer2[k]);
	}

	// Normalise their 3D-parts

	for (int k = 1; k < 4; ++k)
	{
		PauliBuffer1[k] /= std::sqrt(Norm3D_1);
		PauliBuffer2[k] /= std::sqrt(Norm3D_2);
	}

	// Return a final matrix

	return sigma_1 * 0.5*(PauliBuffer1[2]*PauliBuffer2[3] - PauliBuffer1[3]*PauliBuffer2[2]) + 
		   sigma_2 * 0.5*(PauliBuffer1[3]*PauliBuffer2[1] - PauliBuffer1[1]*PauliBuffer2[3]) + 
		   sigma_3 * 0.5*(PauliBuffer1[1]*PauliBuffer2[2] - PauliBuffer1[2]*PauliBuffer2[1]);
}
