#include "Containers.h"
#include <iostream>

int main()
{
	// Testing of Unit

	Unit u1(1, 2, 3, 4);
	Unit u2(2, 3, 4, 5);

	Unit u = u1 * u2;

	std::cout << u(0, 0) << u(0, 1) << u(1, 0) << u(1, 1) << std::endl;

	u(0, 0) *= 5.0;

	std::cout << u(0, 0) << u(0, 1) << u(1, 0) << u(1, 1) << std::endl;

	Unit v1 = u*5.0;
	Unit v2 = 5.0*u;

	std::cout << v1(0, 0) << v1(0, 1) << v1(1, 0) << v1(1, 1) << std::endl;
	std::cout << v2(0, 0) << v2(0, 1) << v2(1, 0) << v2(1, 1) << std::endl;

	Unit v = v1 + v2;

	std::cout << v(0, 0) << v(0, 1) << v(1, 0) << v(1, 1) << std::endl;

	Unit vp = v / 2;

	std::cout << vp(0, 0) << vp(0, 1) << vp(1, 0) << vp(1, 1) << std::endl;

	v /= 2;

	std::cout << v(0, 0) << v(0, 1) << v(1, 0) << v(1, 1) << std::endl;

	Unit* uuu = new Unit(1, 2, 3, 4);

	std::cout << (*uuu)(0, 0) << (*uuu)(0, 1) << (*uuu)(1, 0) << (*uuu)(1, 1) << std::endl;

	Unit sk = skew(sigma_1, sigma_2);

	std::cout << sk(0, 0) << sk(0, 1) << sk(1, 0) << sk(1, 1) << std::endl;	

	// Testing of Vector

	

	return 0;
}