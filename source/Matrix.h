#ifndef __MATRIX_H
#define __MATRIX_H

#include <complex>
typedef std::complex<double> Complex;

class Matrix
// Class which is used for the density matrix and the Hamiltonians
{
public:
	Matrix(Complex e11, Complex e12, Complex e21, Complex e22);
	Matrix(const Matrix& m);
	Matrix(Complex* ptr);
	Matrix();

	// Basic operations over matrices
	Matrix operator+(const Matrix& m) const;
	Matrix operator-(const Matrix& m) const;
	Matrix operator*(const Matrix& m) const;

	// Right multiplication upon a constant
	Matrix operator*(const Complex& Z) const;

	// Left multiplication
	friend Matrix operator*(const Complex& Z, const Matrix& m);
	// 'friend' is used to allow function to see the private elements of the class in 'm'
	// because it's not a member of class

	// Hermite conjugate
	Matrix H() const;

	// Trace function
	Complex Tr() const;

	// Recover the trace = 1 and the hermitance
	void normalise();

	// Double-operations
	void operator+=(const Matrix& m);
	void operator-=(const Matrix& m);
	void operator=(const Matrix& m);
	void operator*=(const Complex& Z); // upon a constant

	// Access method
	Complex get(int i, int j);

	// Change method
	void set(int i, int j, const Complex& Z);

	// Pack/unpack elements to/from an array of 4 elements..
	void pack(Complex* buf);
	void unpack(Complex* buf);

	// .. and only one element from an array
	void put(int i, int j, Complex* buf);

private:
	Complex e11, e12, e21, e22;
};

Matrix comm(const Matrix& m1, const Matrix& m2);
Matrix anticomm(const Matrix& m1, const Matrix& m2);

extern Matrix tau0, tau1, tau2, tau3;

#endif