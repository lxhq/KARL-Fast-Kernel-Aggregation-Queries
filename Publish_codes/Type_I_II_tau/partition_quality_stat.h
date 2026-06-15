#ifndef PARTITION_QUALITY_STAT_H
#define PARTITION_QUALITY_STAT_H

#include "GBF_KC.h"

void sample_Dataset(char*sampleMatrixFileName,double**dataMatrix,double**& sampleMatrix,int dataSize, int sampleSize,int dim);
void offline_KC(char*sampleMatrixFileName,double**dataMatrix,double*alphaArray,int sampleSize,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,char*bulkLoad_TreeName);

#endif