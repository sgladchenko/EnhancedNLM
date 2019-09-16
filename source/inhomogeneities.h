#ifndef __INHOMOGENEITIES_H
#define __INHOMOGENEITIES_H

void generate_random_harmonics(int num_harmonics, double* harmonics);
double mu(double X, double* harmonics, int num_harmonics);

#endif