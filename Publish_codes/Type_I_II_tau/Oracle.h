#ifndef ORACLE_H
#define ORACLE_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "Incremental_Classification.h"

double runOracle(double*q,double**P,svm_node*svm_q,svm_model*model,double*alphaArray,int dim,SVM_stat& stat,svm_node*x_space,int method);
void outputOracle_QueryTime(char*Oracle_TimeFileName,double bestTime);

#endif