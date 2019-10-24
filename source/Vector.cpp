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

Vector Vector::operator+(const Unit& u) const
{

}

Vector Vector::operator-(const Unit& u) const
{

}

Vector Vector::operator*(const Unit& u) const
{


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