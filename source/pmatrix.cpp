#include "pmatrix.h"

PMatrix::PMatrix()
{
	this->e0 = 0; this->e1 = 0; this->e2 = 0; this->e3 = 0;
}

PMatrix::PMatrix(Complex e0, Complex e1, Complex e2, Complex e3)
{
	this->e0 = e0; this->e1 = e1; this->e2 = e2; this->e3 = e3;
}

PMatrix::PMatrix(const PMatrix& pm)
{
	this->e0 = pm.e0; this->e1 = pm.e1; this->e2 = pm.e2; this->e3 = pm.e3;
}

PMatrix::PMatrix(const Matrix& m)
{
	Matrix m1 = tau1 * m;
	Matrix m2 = tau2 * m;
	Matrix m3 = tau3 * m;

	this->e1 = m1.tr();
	this->e2 = m2.tr();
	this->e3 = m3.tr();
	
	this->e0 = m.tr();
}

PMatrix PMatrix::operator+(const PMatrix& pm) const
{
	return PMatrix(this->e0 + pm.e0,
		this->e1 + pm.e1,
		this->e2 + pm.e2,
		this->e3 + pm.e3);
}

PMatrix PMatrix::operator-(const PMatrix& pm) const
{
	return PMatrix(this->e0 - pm.e0,
		this->e1 - pm.e1,
		this->e2 - pm.e2,
		this->e3 - pm.e3);
}

// Right multiplicating on a constant
PMatrix PMatrix::operator*(const Complex& z) const
{
	return PMatrix(this->e0*z, this->e1*z,	this->e2*z, this->e3*z);
}

// Left multiplicating on a constant
PMatrix operator*(const Complex& z, const PMatrix& pm)
{
	return PMatrix(pm.e0*z, pm.e1*z, pm.e2*z, pm.e3*z);
}

void PMatrix::operator=(const PMatrix& pm)
{
	this->e0 = pm.e0;
	this->e1 = pm.e1;
	this->e2 = pm.e2;
	this->e3 = pm.e3;
}

void PMatrix::operator+=(const PMatrix& pm)
{
	this->e0 += pm.e0;
	this->e1 += pm.e1;
	this->e2 += pm.e2;
	this->e3 += pm.e3;
}

void PMatrix::operator-=(const PMatrix& pm)
{
	this->e0 -= pm.e0;
	this->e1 -= pm.e1;
	this->e2 -= pm.e2;
	this->e3 -= pm.e3;
}

void PMatrix::operator*=(const Complex& z)
{
	this->e0 *= z;
	this->e1 *= z;
	this->e2 *= z;
	this->e3 *= z;
}

void PMatrix::set(int i, const Complex& z) // some stupid realisation x2
{
	if (i == 0)
	{
		this->e0 = z;
	}
	else if (i == 1)
	{
		this->e1 = z;
	}
	else if (i == 2)
	{
		this->e2 = z;
	}
	else if (i == 3)
	{
		this->e3 = z;
	}
}

Complex PMatrix::get(int i) // another stupid realisation x2
{
	if (i == 0)
	{
		return this->e0;
	}
	else if (i == 1)
	{
		return this->e1;
	}
	else if (i == 2)
	{
		return this->e2;
	}
	else if (i == 3)
	{
		return this->e3;
	}
}

// Method to form Matrix-objects

Matrix PMatrix::make_matrix() const
{
	return 0.5 * (tau0*this->e0 + tau1*this->e1 + tau2*this->e2 + tau3*this->e3);
}

// Skew-product

PMatrix PMatrix::skew(const PMatrix& pm) const
{
	return PMatrix(0,
		this->e2 * pm.e3 - this->e3 * pm.e2,
		this->e3 * pm.e1 - this->e1 * pm.e3,
		this->e1 * pm.e2 - this->e2 * pm.e1);
}