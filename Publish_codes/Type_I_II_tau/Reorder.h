#ifndef REORDER_H
#define REORDER_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "my_svm.h"

struct orderEntry
{
	int id;
	//method=14
	double dist;
	//method=15
	double norm;
};

//O(nd) time reordering algorithm
void concentric_Reorder(double**P,vector<orderEntry>& orderList,double*alphaArray,int dim,SVM_stat& stat);
void norm_Reorder(double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat);
void oracle_Reorder(double*q,double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat);
void reference_Reorder(double*q,double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat);

void reordering(double**P,double*alphaArray,vector<orderEntry>& orderList,SVM_stat& stat,int dim);

//reordering without creating/ deleting any memory
void init_reorder_mless(double**& PP,double*& reorder_alphaArray,SVM_stat& stat,int dim);
void reordering_mless(double**P,double**PP,double*alphaArray,double*reorder_alphaArray,vector<orderEntry>& orderList,SVM_stat& stat,int dim);
void delete_reorder_mless(double**PP,double*reorder_alphaArray,SVM_stat& stat);

#endif