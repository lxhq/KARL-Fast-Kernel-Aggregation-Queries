#ifndef STAT_ESTIMATION_H
#define STAT_ESTIMATION_H

#include "Reorder.h"
#include "GBF_KC.h"

void stat_Estimation(char*speedUpFileName,int trainNum,double**trainMatrix,double**P,double*alphaArray,SVM_stat& stat,int dim,svm_model*& model, int PT);

void load_sample_idFile(char*sample_idFileName,double**& trainMatrix,int& tNum,double**dataMatrix,int dim);

int loadChunkSize(char*speedUpFileName);
#endif