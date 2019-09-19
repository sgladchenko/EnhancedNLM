#include "matrix.h"

// Pauli matrices
Matrix tau0(1, 0, 0, 1);
Matrix tau1(0, 1, 1, 0);
Matrix tau2(0, -Complex(0, 1.0), Complex(0, 1.0), 0);
Matrix tau3(1, 0, 0, -1);

Matrix::Matrix(Complex e11, Complex e12, Complex e21, Complex e22)
{
	this->e11 = e11; this->e12 = e12; this->e21 = e21; this->e22 = e22;
}

Matrix::Matrix(const Matrix& m)
{
	this->e11 = m.e11;
	this->e12 = m.e12;
	this->e21 = m.e21;
	this->e22 = m.e22;
}

Matrix::Matrix()
{
	this->e11 = 0; this->e12 = 0; this->e21 = 0; this->e22 = 0;
}

Matrix Matrix::operator+(const Matrix& m) const
{
	return Matrix(this->e11 + m.e11,
		this->e12 + m.e12,
		this->e21 + m.e21,
		this->e22 + m.e22
	);
}

Matrix Matrix::operator-(const Matrix& m) const
{
	return Matrix(this->e11 - m.e11,
		this->e12 - m.e12,
		this->e21 - m.e21,
		this->e22 - m.e22
	);
}

Matrix Matrix::operator*(const Matrix& m) const
{
	return Matrix(this->e11 * m.e11 + this->e12 * m.e21,
		this->e11 * m.e12 + this->e12 * m.e22,
		this->e21 * m.e11 + this->e22 * m.e21,
		this->e21 * m.e12 + this->e22 * m.e22
	);
}

// right multiplicating on constant
Matrix Matrix::operator*(const Complex& Z) const
{
	return Matrix(this->e11*Z, this->e12*Z,	this->e21*Z, this->e22*Z);
}

Matrix Matrix::hermit() const
{
	return Matrix(std::conj(this->e11), std::conj(this->e21), std::conj(this->e12), std::conj(this->e22));
}

Complex Matrix::tr() const
{
	return (this->e11 + this->e22);
}

// left multiplicating on constant
Matrix operator*(const Complex& Z, const Matrix& m)
{
	return Matrix(m.e11*Z, m.e12*Z,	m.e21*Z, m.e22*Z);
}

void Matrix::operator=(const Matrix& m)
{
	this->e11 = m.e11;
	this->e12 = m.e12;
	this->e21 = m.e21;
	this->e22 = m.e22;
}

void Matrix::operator+=(const Matrix& m)
{
	this->e11 += m.e11;
	this->e12 += m.e12;
	this->e21 += m.e21;
	this->e22 += m.e22;
}

void Matrix::operator-=(const Matrix& m)
{
	this->e11 -= m.e11;
	this->e12 -= m.e12;
	this->e21 -= m.e21;
	this->e22 -= m.e22;
}

void Matrix::operator*=(const Complex& Z)
{
	this->e11 *= Z;
	this->e12 *= Z;
	this->e21 *= Z;
	this->e22 *= Z;
}

void Matrix::set(int i, int j, const Complex& Z) // some stupid realization
{
	if (i == 0)
	{
		if (j == 0) this->e11 = Z;
		else        this->e12 = Z;
	}
	else
	{
		if (j == 0) this->e21 = Z;
		else        this->e22 = Z;
	}
}

Complex Matrix::get(int i, int j) // another stupid realization
{
	if (i == 0)
	{
		if (j == 0) return this->e11;
		else        return this->e12;
	}
	else
	{
		if (j == 0) return this->e21;
		else        return this->e22;
	}
}

Matrix comm(const Matrix& m1, const Matrix& m2)
{
	return m1*m2 - m2*m1;
}

Matrix anticomm(const Matrix& m1, const Matrix& m2)
{
	return m1*m2 + m2*m1;
}

void Matrix::pack(Complex* buf)
{
	buf[0] = this->e11;
	buf[1] = this->e12;
	buf[2] = this->e21;
	buf[3] = this->e22;
}

void Matrix::unpack(Complex* buf)
{
	this->e11 = buf[0];
	this->e12 = buf[1];
	this->e21 = buf[2];
	this->e22 = buf[3];
}

void Matrix::put(int i, int j, Complex* buf)
{
	if (i == 0)
	{
		if (j == 0) buf[0] = this->e11;
		else        buf[0] = this->e12;
	}
	else
	{
		if (j == 0) buf[0] = this->e21;
		else        buf[0] = this->e22;
	}
}

void Matrix::normalize()
{
	// recover hermitance
	Complex o11, o12, o21, o22; // old ones

	o11 = this->e11;
	o12 = this->e12;
	o21 = this->e21;
	o22 = this->e22;

	this->e11 = 0.5*(o11 + std::conj(o11));
	this->e12 = 0.5*(o12 + std::conj(o21));
	this->e21 = 0.5*(o21 + std::conj(o12));
	this->e22 = 0.5*(o22 + std::conj(o22));

	// recover trace = 1
	Complex tr = this->e11 + this->e22;

	this->e11 /= tr;
	this->e12 /= tr;
	this->e21 /= tr;
	this->e22 /= tr;
}