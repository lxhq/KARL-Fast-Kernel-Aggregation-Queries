#ifndef MY_SVM_H
#define MY_SVM_H

#include "init_KC.h"
#include "svm.h"

void convert_2D_qMatrix_to_SVMFormat(double**queryMatrix,svm_node**& svm_qMatrix,int dim,int qNum);
void loadSVM_model(char*modelName,svm_model*& model);
void prediction(const svm_model *model, const svm_node *x,SVM_stat& stat);
void reorder_model(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model);
void init_model_inMemory(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space);
void createModel_inMemory(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space);

#endif