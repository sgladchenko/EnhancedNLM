#include "Containers.h"

Layer::Layer()
{
	MyCountX = 0;
}

Layer::Layer(int given_MyCountX, Complex* buffer)
{
	MyCountX = given_MyCountX;
	VectorBuffer = new Vector[MyCountX+2]; // MyCountX+2: to consider two additional points on the left and on the right

	for (int x = 0; x < MyCountX+2; ++x)
	{
		VectorBuffer[x].Init(buffer + x*16*N_E);
	}
}

void Layer::Init(int given_MyCountX, Complex* buffer)
{
	// If it was not empty
	if (MyCountX != 0)
	{
		delete [] VectorBuffer;
	}

	MyCountX = given_MyCountX;
	VectorBuffer = new Vector[MyCountX+2]; // MyCountX+2: to consider two additional points on the left and on the right

	for (int x = 0; x < MyCountX+2; ++x)
	{
		VectorBuffer[x].Init(buffer + x*16*N_E);
	}
}

Vector& Layer::operator()(int x)
{
	return VectorBuffer[x];
}

void Layer::Dump(Complex* buf)
{
	for (int x = 1; x <= MyCountX; ++x) // Dump only its own part of the whole grid of x-coordinates
	{
		VectorBuffer[x].Pack(buf + (x-1)*16*N_E);
	}
}

Layer::~Layer()
{
	delete [] VectorBuffer;
}