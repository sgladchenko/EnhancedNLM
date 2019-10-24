#include "Log.h"

Log::Log(int given_MyRank, double given_Start_time, double given_Z_init, double given_Z_displacement)
{
	_Start_time = given_Start_time;
	_MyRank = given_MyRank;

	_Z_init = given_Z_init;
	_Z_displacement = given_Z_displacement;
}

void Log::OutTime(std::string message)
{
	if (_MyRank == 0)
	{
		std::cout << message << "; elapsed=" << MPI_Wtime() - _Start_time << " (sec)" << std::endl;
	}
}

void Log::Out(std::string message)
{
	if (_MyRank == 0)
	{
		std::cout << message << std::endl;
	}
}

void Log::OutIter(double L)
{
	if (_MyRank == 0)
	{
		std::cout << std::fixed << std::setprecision(2) << "z=" << L/ph::km << " of "
		          << (_Z_init + _Z_displacement)/ph::km << " (km); elapsed="
		          << MPI_Wtime() - _Start_time << " (sec)" << std::endl;
	}
}

void Log::OutDebug(double rank, double smth)
{
	if (_MyRank == rank)
	{
		std::cout << smth << std::endl;
	}
}
