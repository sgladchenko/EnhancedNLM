#ifndef __STANDARDK_H
#define __STANDARDK_H

#include "Matrix.h"
#include "Constants.h"

void Interchange_Ks(Matrix*  Kbuf,
					Complex* LeftSend_buf,
					Complex* RightSend_buf,
					Complex* LeftRecv_buf,
					Complex* RightRecv_buf,
					int      LeftNeighbour,
					int      RightNeighbour,
					int      Korder,
					int      MyCountX);

void Evaluate_K1(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K1); // the output buffer

void Evaluate_K2(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K1,
				 Matrix* K2); // the output buffer

void Evaluate_K3(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K2,
				 Matrix* K3); // the output buffer

void Evaluate_K4(Matrix* old_p,
				 Matrix* old_m,
				 Matrix* old_antip,
				 Matrix* old_antim,
				 Matrix* H_0_mesh,
				 double* Spec_nu,
				 double* Spec_antinu,
				 double* Harmonics,
				 int     MyCountX,
				 double  Z,
				 double  Xleft,
				 Matrix* K3,
				 Matrix* K4); // the output buffer

#endif