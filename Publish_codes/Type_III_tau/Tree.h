#ifndef TREE_H
#define TREE_H

#include "init_KC.h"
#include "Euclid_Bound.h"
#include "my_svm.h"

class Node;
class Tree;

class Node
{
public:
	double sum_alpha_pos;
	double sum_alpha_neg;
	vector<int> idList; //idxs
	vector<Node*> childVector;
	svm_model* model;
	svm_node* x_space;

	virtual double LB(double*q,int dim,SVM_stat& stat)=0;
	virtual double UB(double*q,int dim,SVM_stat& stat)=0;

	virtual Node*createNode()=0;
	virtual void update_Aug(Node*node,Tree*t)=0;

	//used in leaf node
	void initModel_inMemory(double**dataMatrix,double*outputArray,int dim,SVM_stat& stat);
	void createModel_inMemory(double**dataMatrix,double*outputArray,int dim,SVM_stat& stat);
	void freeModel_inMemory();
};

class Tree
{
public:
	int dim;
	double**dataMatrix;
	double*outputArray;
	int leafCapacity; //Set to be 20 in tKDC
	SVM_stat stat;
	Node*rootNode;
};

#endif