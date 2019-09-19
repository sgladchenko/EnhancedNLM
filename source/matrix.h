#ifndef __MATRIX_H
#define __MATRIX_H

#include <complex>
typedef std::complex<double> Complex;

// class which is used for density matrix
class Matrix
{
public:
	Matrix(Complex e11, Complex e12, Complex e21, Complex e22);
	Matrix(const Matrix& m);
	Matrix();

	Matrix operator+(const Matrix& m) const;
	Matrix operator-(const Matrix& m) const;
	Matrix operator*(const Matrix& m) const;
	Matrix operator*(const Complex& Z) const;

	// Hermitian conjugate
	Matrix hermit() const;

	// trace function
	Complex tr() const;

	void normalize(); // hermitance + trace=1

	void operator+=(const Matrix& m);
	void operator-=(const Matrix& m);
	void operator*=(const Complex& Z);
	void operator=(const Matrix& m);

	friend Matrix operator*(const Complex& Z, const Matrix& m); // left multiplicating
	// friend is used to allow function to see the private elements of the class in 'm'
	// because it's not a member of class

	void set(int i, int j, const Complex& Z);
	Complex get(int i, int j);

	// pack elements to an array of 4 elements
	void pack(Complex* buf);
	void unpack(Complex* buf);

	// and only one element
	void put(int i, int j, Complex* buf);

private:
	Complex e11, e12, e21, e22;
};

Matrix comm(const Matrix& m1, const Matrix& m2);
Matrix anticomm(const Matrix& m1, const Matrix& m2);

extern Matrix tau0, tau1, tau2, tau3;

#endif