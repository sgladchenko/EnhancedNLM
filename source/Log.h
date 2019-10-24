#ifndef __LOGH
#define __LOGH

#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <string>

#include "Constants.h"

class Log
{
	public:
		Log(int given_MyRank, double given_Start_time, double given_Z_init, double given_Z_displacement);

		// Different ways of output
		void OutTime(std::string message);
		void Out(std::string message); // simply a message without a new time stamp
		void OutIter(double L);

		// I wish I'd never use it...
		void OutDebug(double rank, double smth);

	protected:
		int _MyRank;
		double _Start_time;

		// These values will be shown with the iteration stamps
		double _Z_init;
		double _Z_displacement;
};

#endif