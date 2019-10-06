#ifndef __SEGMENTATION_H
#define __SEGMENTATION_H

void get_segment_X(int rank, int size, int& left, int& right);
void get_countXs(int size, int* countXs);

void get_scatter_counts(int size, int* counts);
void get_scatter_displacements(int size, int* displacements);

void get_line_counts(int size, int* counts);
void get_line_displacements(int size, int* displacements);

void get_ad_counts(int size, int* ad_counts);
void get_ad_displacements(int size, int* ad_displacements);

#endif