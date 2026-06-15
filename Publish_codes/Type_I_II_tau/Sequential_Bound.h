#ifndef SEQUENTIAL_BOUND_H
#define SEQUENTIAL_BOUND_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "SS.h"
#include "my_svm.h"

//method=9,11
void pre_Compute_sequential(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double*& a_G,double& S_G,double*& center,double& radius,SVM_stat& stat);

//method=9/10
void sequential_bound(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,int method);

//method=30/31
void pre_Compute_sequential_MBR(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double**& boundary,SVM_stat& stat);
void sequential_bound_MBR(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double**boundary,SVM_stat& stat,int method);

//method=32/33
void pre_Compute_sequential_Delta(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double*& center,double& radius,SVM_stat& stat);
void sequential_bound_Delta(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double*center,double radius,SVM_stat& stat,int method);

//method=11
void sequential_bound_sparse(double*q,double*alphaArray,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,double*p_sqNormArray,double**lu_Table);

#endif