#ifndef INCREMENTAL_CLASSIFICATION_H
#define INCREMENTAL_CLASSIFICATION_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "my_svm.h"
#include "Reorder.h"

struct group_Info
{
	double w_G;
	double*a_G;
	double S_G;
	double*center;
	double radius;
};

double LB_T(double*q,group_Info& g_Info,double& gamma_sum,int dim,SVM_stat& stat);
double UB_C(double*q,group_Info& g_Info,double& gamma_sum,int dim,SVM_stat& stat);
void incremental_Classification(double*q,double**P,group_Info*groupArray,double*alphaArray,int chunkSize,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q);
void incremental_Classification_multiStep(double*q,double**P,group_Info*groupArray,double*alphaArray,int chunkSize,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q);

void create_groupArray(double**P,double*alphaArray,group_Info*& groupArray,int chunkSize,int dim,SVM_stat& stat,bool initPruneCounter);
#endif