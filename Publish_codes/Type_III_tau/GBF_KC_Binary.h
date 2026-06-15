#ifndef GBF_KC_BINARY_H
#define GBF_KC_BINARY_H

#include "init_KC.h"
#include "SS_Binary.h"
#include "kd_tree.h"
#include "m_tree.h"
#include "my_svm.h"

struct pqNode
{
	Node*node;
	double discrepancy;
	double node_L;
	double node_U;
};

//This is the maximum heap
struct comparePriority 
{
	bool operator()(pqNode& p1, pqNode& p2)
	{
		return p1.discrepancy<p2.discrepancy;
	}
};

typedef priority_queue<pqNode,vector<pqNode>,comparePriority> PQ;

double computeSqNorm(double*q,int dim);
void GBF_iter(double*q,svm_node*svm_q,Tree& tree,int dim,SVM_stat& stat);
void KC_Algorithm(double**queryMatrix,double**dataMatrix,double*outputArray,int qNum,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,char*bulkLoad_TreeName);

#endif