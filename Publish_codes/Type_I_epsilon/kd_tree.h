#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"

class kdNode : public Node
{
public:
	double**boundary; //boundary

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);

	void update_Aug(Node*node,Tree*t){}

	kdNode*createNode();
};

//kd-tree +augment Tri-Inequality Information
class kdAugNode : public kdNode
{
public:
	//Augment Info
	double*center;
	double radius;

	//facilitates online (sharing) computation
	double temp_obt_dist;

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);

	void updateAugInfo(kdAugNode*node,Tree*t);
	void update_Aug(Node*node,Tree*t);

	kdAugNode*createNode();
};

class kdTree : public Tree
{
public:
	//Member functions
	kdTree(int dim,double**dataMatrix,int leafCapacity,KDE_stat& stat);

	void getNode_Boundary(kdNode*node);
	double obtain_SplitValue(kdNode*node,int split_Dim);
	void KD_Tree_Recur(kdNode*node,int split_Dim);
	void build_kdTree(KDE_stat& stat);
	void initTree_sumE(kdNode*node);

	void updateAugment(kdNode*node);

	//added this function for testing on 27/4/2018 
	int select_split_Dim(kdNode*node);
};

class kdLinearAugNode : public kdAugNode
{
public:
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);

	void update_linearAugInfo(kdLinearAugNode*node,Tree*t);

	void update_Aug(Node*node,Tree*t);
	kdLinearAugNode*createNode();
};

#endif