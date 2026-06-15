#ifndef CACHE_ORDER_H
#define CACHE_ORDER_H

#include "init_KC.h"
#include "Reorder.h"
#include "binaryTree.h"

struct pivotSet
{
	double**pivotMatrix;
	int num_of_pivots;
};

//offline
void init_pivotSet(char*pivotFileName,pivotSet& pSet,int dim);
void build_multiple_Tree(vector<binaryTree>& multi_b_Tree,pivotSet& pSet,int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat,int method);

//online
int NN(double*q,pivotSet& pSet,int dim);

#endif