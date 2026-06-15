#ifndef WEAK_AUTO_H
#define WEAK_AUTO_H

#include "Sequential_Bound.h"
#include "GBF_KC.h"

int weak_auto(double**queryMatrix,double**dataMatrix,double*alphaArray,svm_node**svm_qMatrix,svm_model*& model,int sNum,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,int leafCapacity,kdTree_adv& kd_Tree_adv);

#endif