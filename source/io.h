#ifndef __IOFUNCTIONS_H
#define __IOFUNCTIONS_H

#include "constants.h"
#include <fstream>

// Filenames of binary output

#define P_OUTBIN "p.bin"
#define M_OUTBIN "m.bin"
#define ANTIP_OUTBIN "antip.bin"
#define ANTIM_OUTBIN "antim.bin"

#define REC_FNAME "rec.bin"

#define ZGRID "gridz.bin"
#define XGRID "gridx.bin"
#define EGRID "gridE.bin"

#define Z_FNAME "z.txt"

#define HARMONICS_FILE "harmonics.bin"

void Text_input_Z_init(double& Z_init);
void Text_output_Z_final(double Z_final);

void Bin_input_scatt_buffer(int size, Complex* scatt_buffer);
void Bin_output_rec(int size, int* counts, int* displacements, Complex* gath_buffer);
void Bin_output_line(int size, int* counts, int* displacements, Complex* gath_buffer);
void Bin_input_harmonics(int num_harmonics, double* harmonics);

void make_gridx();
void make_gridE();
void make_gridz();
void make_bin_files();
void add_to_gridz(double Z_here);

#endif