#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"

class kdNode : public Node
{
public:
	double**boundary_pos; //boundary_pos
	double**boundary_neg; //boundary_neg

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	void update_Aug(Node*node,Tree*t){}

	kdNode*createNode();
};

class kdLinearAugNode : public kdNode
{
public:
	double*a_G_pos;
	double*a_G_neg;
	double S_G_pos;
	double S_G_neg;

	//facilitates online (sharing) computation
	double gamma_sum_pos;
	double gamma_sum_neg;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	void update_linearAugInfo(kdLinearAugNode*node,Tree*t);

	void update_Aug(Node*node,Tree*t);
	kdLinearAugNode*createNode();
};

class kdLinearAugNode_adv : public kdLinearAugNode
{
public:
	kdLinearAugNode_adv*createNode();
	void update_linearAugInfo_adv(kdLinearAugNode*node,Tree*t);
	void update_Aug(Node*node,Tree*t);
};

class kdTree : public Tree
{
public:
	//Member functions
	kdTree(int dim,double**dataMatrix,double*outputArray,int leafCapacity,SVM_stat& stat);

	void getNode_Boundary(kdNode*node);
	double obtain_SplitValue(kdNode*node,int split_Dim);
	void KD_Tree_Recur(kdNode*node,int split_Dim);
	void build_kdTree(SVM_stat& stat);
	void initTree_alpha(kdNode*node);

	void updateAugment(kdNode*node);
};

class kdTree_adv : public kdTree
{
public:
	//Member functions
	kdTree_adv(int dim,double**dataMatrix,double*outputArray,int leafCapacity,SVM_stat& stat)
		: kdTree(dim,dataMatrix,outputArray,leafCapacity,stat){}
	void KD_Tree_adv_Recur(kdNode*node,int split_Dim);
	void build_kdTree_adv(SVM_stat& stat);
	void merge_Boundary(kdNode*node,double**boundary_1,double**boundary_2,bool isPos);
};

#endif