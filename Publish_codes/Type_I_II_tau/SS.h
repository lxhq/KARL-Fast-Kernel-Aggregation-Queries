#ifndef SS_H
#define SS_H

#include "init_KC.h"
#include "inverted_Index.h"
#include "exp_LookUp.h"

//handle sparse vector via inverted index (method=11)
void SS_sparse_iter(double*q,double*alphaArray,int dim,SVM_stat& stat,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,double q_sqNorm,double*p_sqNormArray,double**lu_Table);

//handle dense vector
void SS_iter(double*q,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);

//This is used for testing the sequential scanning
void convert_to_1d(double**dataMatrix,double*&dataMatrix_1d,int dim,SVM_stat& stat);
void SS_iter_1d(double*q,double*dataMatrix_1d,double*alphaArray,int dim,SVM_stat& stat);

#endif