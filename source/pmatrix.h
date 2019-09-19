#ifndef __PMATRIX_H
#define __PMATRIX_H

#include "matrix.h"

class PMatrix // Implementation of the Matrix, but in terms of Pauli-components and for U(2)
{
public:
	PMatrix();
	PMatrix(Complex e0, Complex e1, Complex e2, Complex e3); // matrix = e_a t^a; t^a = {1/2, tau^k / 2}
	PMatrix(const PMatrix& pm);
	PMatrix(const Matrix& m);

	// Note: in futher e_k are presumed to be real-valued; Complex has been put here to avoid
	//       an incovenience, caused by 'rough' calculations

	PMatrix operator+(const PMatrix& pm) const;
	PMatrix operator-(const PMatrix& pm) const;
	PMatrix operator*(const Complex& z) const;

	// Pull out the corresponding Matrix object
	Matrix make_matrix() const;

	// Skew-symmetric product
	PMatrix skew(const PMatrix& pm) const; // e_{abc} this^{b} pm^{c}; a,b,c > 0 (w/out the unit)

	friend PMatrix operator*(const Complex& z, const PMatrix& pm); // left multiplicating
	// friend is used to allow function to see private elements of the class in 'm'
	// because it's not a member of class

	// Some other operations
	void operator=(const PMatrix& pm);
	void operator+=(const PMatrix& pm);
	void operator-=(const PMatrix& pm);
	void operator*=(const Complex& z);

	// Access methods
	void set(int i, const Complex& z);
	Complex get(int i);

private:
	Complex e0, e1, e2, e3;
};

#endif