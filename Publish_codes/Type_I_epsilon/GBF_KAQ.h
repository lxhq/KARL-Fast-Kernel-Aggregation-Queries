#ifndef GBF_KC_H
#define GBF_KC_H

#include "init_KAQ.h"
#include "SS.h"
#include "kd_tree.h"
#include "m_tree.h"
#include "binaryTree.h"
#include "Validation.h"

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

void GBF_iter(double*q,Tree& tree,int dim,KDE_stat& stat);

void KAQ_Algorithm(double**queryMatrix,double**dataMatrix,int qNum,int dim,int leafCapacity,int internalCapacity,int method,KDE_stat& stat,char*bulkLoad_TreeName);

#endif