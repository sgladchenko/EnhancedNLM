#include "Constants.h"

class Unit
// Actually it's a re-issue of the Matrix-class -- rho(type, zeta, x_k, z_n, E)
{
	public:
		Unit(Complex* buffer);
		Unit(Complex e0, Complex e1, Complex e2, Complex e3);
		Unit(const Unit& u);
		Unit();

		void Init(Complex e0, Complex e1, Complex e2, Complex e3);
		void Init(Complex* buffer);

		// Operations over these objects

		// Singular ones:
		Unit operator+(const Unit& u) const;
		Unit operator-(const Unit& u) const;
		Unit operator*(const Unit& u) const;

		Unit operator/(const Complex& z) const;

		// Right multiplication on a Complex number
		Unit operator*(const Complex z) const;

		// Left multiplication on a Complex number
		friend Unit operator*(const Complex& z, const Unit& u);

		// Double-operations
		void operator+=(const Unit& u);
		void operator-=(const Unit& u);
		void operator*=(const Complex& z);
		void operator/=(const Complex& z);

		void operator=(const Unit& u);

		// Other matrix stuff

		// Trace function
		Complex Tr() const;

		// Recover trace = 1 and Hermitance
		void normalise();

		// Return the norm of the matrix
		Complex norm();

		// Access/modify methods
		Complex& operator()(int i, int j);
		Complex& operator[](int k);

		// Commutator and the skew-product
		friend Unit commutator(const Unit& u1, const Unit& u2);
		friend Unit skew(const Unit& u1, const Unit& u2);

	protected:
		Complex ComplexBuffer[4];
};

// Pauli matrices declaration
extern Unit Id, sigma_1, sigma_2, sigma_3;
extern Unit PauliMatrices[4];

class Vector
// The implementation of an object of the type rho(x_k, z_n)
{
	public:
		Vector(Complex* buffer);
		Vector();

		void Init(Complex* buffer);

		// Singular operations
		Vector operator+(const Vector& v) const;
		Vector operator-(const Vector& v) const;
		Vector operator*(const Vector& v) const;

		Vector operator/(const Complex& z) const;

		// Right multiplication on a Complex number
		Vector operator*(const Complex z) const;

		// Left multiplication on a Complex number
		friend Vector operator*(const Complex& z, const Vector& v);

		// Double-operations
		void operator+=(const Vector& v);
		void operator-=(const Vector& v);
		void operator*=(const Complex& z);
		void operator/=(const Complex& z);

		void operator=(const Vector& v);

		// Access/modify method
		Unit& operator()(int e, int type, int zeta);

		// Calculational stuff

		Complex norm();

		// Recover trace = 1 and Hermitance of each matrix
		void normalise();

		// Commutator and the skew-product
		friend Vector commutator(const Vector& v1, const Vector& v2);
		friend Vector skew(const Vector& v1, const Vector& v2);

	protected:
		Unit UnitBuffer[N_E*4];
};

class Layer
// The implementation of an object to be operated as a layer of matrices rho(z_n)
{
	public:
		// Initialise from a given ScatterBuffer pointer
		Layer(int given_MyCountX, Complex* buffer);
		Layer();

		void Init(int given_MyCountX, Complex* buffer);

		Vector& operator()(int x);

		~Layer();

	protected:
		// Number of points in this chunck
		int MyCountX;

		// Main buffers of Complex-valued numbers, containing elements of the matrices
		Vector* VectorBuffer;
};