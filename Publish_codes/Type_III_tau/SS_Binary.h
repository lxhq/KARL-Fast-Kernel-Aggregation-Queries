#ifndef SS_BINARY_H
#define SS_BINARY_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "my_svm.h"

//Without LibSVM
void SS_iter(double*q,double**dataMatrix,double*outputArray,int dim,SVM_stat& stat);

//LibSVM_{B+}
void pre_Compute_Bplus(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double*& a_G_pos,double*& a_G_neg,double& S_G_pos,double& S_G_neg,double*& center_pos,double*& center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat);
void LibSVM_Bplus(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double*& a_G_pos,double*& a_G_neg,double& S_G_pos,double& S_G_neg,double*center_pos,double*center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat);

//LibSVM_{MBR}
void pre_Compute_MBR(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double**& boundary_pos,double**& boundary_neg,SVM_stat& stat);
void LibSVM_MBR(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double**& boundary_pos,double**& boundary_neg,SVM_stat& stat);
//LibSVM_{Delta}
void pre_Compute_Delta(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double*& center_pos,double*& center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat);
void LibSVM_Delta(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double*center_pos,double*center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat);
#endif