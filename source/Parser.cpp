#include <fstream>
#include <string>

#include "BaseBin.h"
#include "Constants.h"

// Some templates to simplify writing/reading binaries

template<typename T>
void bin_read(std::ifstream& stream, T* buf, int count)
{
	stream.read(reinterpret_cast<char*>(buf), count*sizeof(T));
}

template<typename T>
void bin_seek(std::ifstream& stream, int count)
{
	stream.seekg(count*sizeof(T), std::ios::beg);
}

template<typename T>
void bin_write(std::ofstream& stream, T* buf, int count)
{
	stream.write(reinterpret_cast<char*>(buf), count*sizeof(T));
}

// Get parameters of grids

void GetGridNumbers(int& localN_Z, int& localN_X, int& localN_E)
{
	std::ifstream gridz, gridx, gridE;

	gridz.open(ZGRID, std::ifstream::binary);
	gridx.open(XGRID, std::ifstream::binary);
	gridE.open(EGRID, std::ifstream::binary);

	gridz.seekg(0, gridz.end);
	gridx.seekg(0, gridx.end);
	gridE.seekg(0, gridE.end);

	localN_Z = gridz.tellg() / sizeof(double);
	localN_X = gridx.tellg() / sizeof(double);
	localN_E = gridE.tellg() / sizeof(double);

	gridz.close();
	gridx.close();
	gridE.close();
}

// 2D output //

void Fix_E(std::string fin, std::string fout, int e) // e is the number of bin (from 0 to N_E)
{
	std::ifstream fullbin;
	fullbin.open(fin, std::ifstream::binary);

	std::ofstream outbin;
	outbin.open(fout, std::ofstream::binary);

	int localN_Z, localN_X, localN_E;
	GetGridNumbers(localN_Z, localN_X, localN_E);

	Complex tmp[4];

	for (int z = 0; z < localN_Z; ++z)
	{
		for (int x = 0; x < localN_X; ++x)
		{
			bin_seek<Complex>(fullbin, z*localN_X*localN_E*4 + x*localN_E*4 + e*4);
			bin_read(fullbin, tmp, 4);
			bin_write(outbin, tmp, 4);
		}
	}

	fullbin.close();
	outbin.close();
}

// 1D output //

void Average_by_x(std::string fin, std::string fout)
{
	std::ifstream fullbin;
	fullbin.open(fin, std::ifstream::binary);

	std::ofstream outbin;
	outbin.open(fout, std::ofstream::binary);

	int localN_Z, localN_X, localN_E;
	GetGridNumbers(localN_Z, localN_X, localN_E);

	Complex* tmp = new Complex[4*localN_E];
	Complex* buf = new Complex[4*localN_E];

	for (int z = 0; z < localN_Z; ++z)
	{
		for (int i = 0; i < 4*localN_E; ++i)
		{
			tmp[i] = 0;
		}

		for (int x = 0; x < localN_X; ++x)
		{
			bin_seek<Complex>(fullbin, z*localN_X*localN_E*4 + x*localN_E*4);
			bin_read(fullbin, buf, 4*localN_E);

			for (int e = 0; e < localN_E; ++e)
			{
				for (int i = 0; i < 4; ++i)
				{
					tmp[4*e + i] += buf[4*e + i];
				}
			}
		}

		for (int i = 0; i < 4*localN_E; ++i)
		{
			tmp[i] /= localN_X;
		}		

		bin_write(outbin, tmp, 4*localN_E);
	}

	fullbin.close();
	outbin.close();
}

int main()
{
	Fix_E(P_BIN, "./data/Fix_E.bin", 3);
	Average_by_x(P_BIN, "./data/Average_by_x.bin");
}