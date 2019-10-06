#ifndef __IOFUNCTIONS_H
#define __IOFUNCTIONS_H

#include "constants.h"
#include <fstream>

// Filenames of binary output

#define P_OUTBIN "./data/p.bin"
#define M_OUTBIN "./data/m.bin"
#define ANTIP_OUTBIN "./data/antip.bin"
#define ANTIM_OUTBIN "./data/antim.bin"

#define REC_FNAME "./data/rec.bin"

#define ZGRID "./data/gridz.bin"
#define XGRID "./data/gridx.bin"
#define EGRID "./data/gridE.bin"

#define Z_FNAME "./data/z.txt"

#define HARMONICS_FILE "./data/harmonics.bin"

#define AD_OUT "./data/ad.bin"

void Text_input_Z_init(double& Z_init);
void Text_output_Z_final(double Z_final);

void Bin_input_scatt_buffer(int size, Complex* scatt_buffer);
void Bin_output_rec(int size, int* counts, int* displacements, Complex* gath_buffer);
void Bin_output_line(int size, int* counts, int* displacements, Complex* gath_buffer);
void Bin_input_harmonics(int num_harmonics, double* harmonics);

void Bin_output_ad(int size, int* ad_counts, int* ad_displacements, Complex* adiabaticity_buffer);

void make_gridx();
void make_gridE();
void make_gridz();
void make_bin_files();
void add_to_gridz(double Z_here);

#endif