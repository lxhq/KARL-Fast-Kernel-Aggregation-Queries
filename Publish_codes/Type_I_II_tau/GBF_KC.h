#ifndef GBF_KC_H
#define GBF_KC_H

#include "init_KC.h"
#include "SS.h"
#include "kd_tree.h"
#include "m_tree.h"
#include "binaryTree.h"
#include "Sequential_Bound.h"
#include "my_svm.h"
#include "Incremental_Classification.h"
#include "Stat_Estimation.h"
#include "Oracle.h"
#include "Oracle_Tree.h"
#include "cache_order.h"
//#include "weak_auto.h"

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

//The tree can be any tree
void GBF_iter(double*q,svm_node*svm_q,Tree& tree,int dim,SVM_stat& stat);

double computeSqNorm(double*q,int dim);

void KC_Algorithm(double**queryMatrix,double**dataMatrix,double*alphaArray,int qNum,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,svm_model*& model,int chunkSize,char*& oracle_TimeFileName,char*pivotFileName,char*bulkLoad_TreeName);

#endif