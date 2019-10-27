#include "Containers.h"

// Initialisers

Vector::Vector(Complex* buffer)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta].Init(buffer + 4*(e*4 + type*2 + zeta));
			}
		}
	}
}

Vector::Vector()
{
	for (int k = 0; k < 4*N_E; ++k)
	{
		UnitBuffer[k].Init(0, 0, 0, 0);
	}
}

void Vector::Init(Complex* buffer)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta].Init(buffer + 4*(e*4 + type*2 + zeta));
			}
		}
	}
}

// Singular operations

Vector Vector::operator+(const Vector& v) const
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = UnitBuffer[e*4 + type*2 + zeta] + v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}

	return result;
}

Vector Vector::operator-(const Vector& v) const
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = UnitBuffer[e*4 + type*2 + zeta] - v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}

	return result;
}

Vector Vector::operator*(const Vector& v) const
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = UnitBuffer[e*4 + type*2 + zeta] * v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}

	return result;
}

Vector Vector::operator/(const Complex& z) const
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = UnitBuffer[e*4 + type*2 + zeta] / z;
			}
		}
	}

	return result;	
}

// Right multiplication on a Complex number
Vector Vector::operator*(const Complex z) const
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = UnitBuffer[e*4 + type*2 + zeta] * z;
			}
		}
	}

	return result;	
}

// Left multiplication on a Complex number
Vector operator*(const Complex& z, const Vector& v)
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = v.UnitBuffer[e*4 + type*2 + zeta] * z;
			}
		}
	}

	return result;
}

// Double-operations

void Vector::operator+=(const Vector& v)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta] += v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}
}

void Vector::operator-=(const Vector& v)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta] -= v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}
}

void Vector::operator*=(const Complex& z)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta] *= z;
			}
		}
	}
}

void Vector::operator/=(const Complex& z)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta] /= z;
			}
		}
	}
}

void Vector::operator=(const Vector& v)
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta] = v.UnitBuffer[e*4 + type*2 + zeta];
			}
		}
	}
}

// Access & modify method
Unit& Vector::operator()(int e, int type, int zeta)
{
	return UnitBuffer[e*4 + type*2 + zeta];
}

// Other calculational stuff

Complex Vector::norm()
{
	Complex Norm = 0;

	for (int k = 0; k < 4*N_E; ++k)
	{
		Norm += UnitBuffer[k].norm();
	}

	return Norm;
}

// Recover trace = 1 and Hermitance of each matrix
void Vector::normalise()
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				UnitBuffer[e*4 + type*2 + zeta].normalise();
			}
		}
	}
}

// Commutator and the skew-product
Vector commutator(const Vector& v1, const Vector& v2)
{
	return v1*v2 - v2*v1;
}

Vector skew(const Vector& v1, const Vector& v2)
{
	Vector result;

	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			for (int zeta = 0; zeta < 2; ++zeta)
			{
				result(e, type, zeta) = skew(v1.UnitBuffer[e*4 + type*2 + zeta], v2.UnitBuffer[e*4 + type*2 + zeta]);
			}
		}
	}

	return result;
}

// Multiply components with different zeta by \pm
void Vector::MultiplyByPrefactor()
{
	for (int e = 0; e < N_E; ++e)
	{
		for (int type = 0; type < 2; ++type)
		{
			// For zeta = 1 (this's equal to $\zeta=-1$) we should multiply by -1, if needed
			UnitBuffer[e*4 + type*2 + 1] *= -1.0;
		}
	}	
}

// Pack/unpack from the external Complex-buffer
void Vector::Pack(Complex* buf)
{
	for (int k = 0; k < 4*N_E; ++k)
	{
		UnitBuffer[k].Pack(buf + k*4);
	}
}

void Vector::Unpack(Complex* buf)
{
	for (int k = 0; k < 4*N_E; ++k)
	{
		UnitBuffer[k].Unpack(buf + k*4);
	}
}