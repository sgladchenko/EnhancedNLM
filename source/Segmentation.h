#ifndef __SEGMENTATION_H
#define __SEGMENTATION_H

void Obtain_SegmentX(int MyRank, int CommSize, int& MyLeftX, int& MyRightX);
void Obtain_CountXs(int CommSize, int* CountXs);

void Obtain_ScatterCounts(int CommSize, int* ScatterCounts);
void Obtain_ScatterDisplacements(int CommSize, int* ScatterDisplacements);

void Obtain_GatherCounts(int CommSize, int* GatherCounts);
void Obtain_GatherDisplacements(int CommSize, int* GatherDisplacements);

#endif