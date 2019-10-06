#ifndef __EVOLUTION_H
#define __EVOLUTION_H

#include "matrix.h"
#include "constants.h"

void interchange_Ks(Matrix*  Kbuf,
					Complex* leftsend,
					Complex* rightsend,
					Complex* leftrecv,
					Complex* rightrecv,
					int left_neighbour,
					int right_neighbour,
					int Korder,
					int mycountX);

void interchange_PH_p(PMatrix* PH_p,
					  Complex* send_PH_p_buf,
					  Complex* recv_PH_p_buf,
					  int left_neighbour,
					  int right_neighbour,
					  int mycountX);

void evaluate_adiabaticity(PMatrix* PH_p_old,
						   PMatrix* PH_p,
						   Complex* adiabaticity,
						   int mycountX);

void evaluate_K1(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 PMatrix* PH_p,
				 Matrix* K1); // the output buffer

void evaluate_K2(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K1,
				 Matrix* K2); // the output buffer

void evaluate_K3(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K2,
				 Matrix* K3); // the output buffer

void evaluate_K4(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* harmonics,
				 int mycountX,
				 double Z,
				 double Xleft,
				 Matrix* K3,
				 Matrix* K4); // the output buffer

#endif