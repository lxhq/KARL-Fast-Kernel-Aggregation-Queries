#ifndef KD_TREE_H
#define KD_TREE_H

#include "Tree.h"

class kdNode : public Node
{
public:
	double**boundary; //boundary

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

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

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	void updateAugInfo(kdAugNode*node,Tree*t);
	void update_Aug(Node*node,Tree*t);

	kdAugNode*createNode();
};

class kdAugNode_MIN : public kdNode
{
public:
	//Augment Info
	double**d_SortObj;
	double UB(double*q,int dim,SVM_stat& stat);

	void updateAugMinInfo(kdAugNode_MIN*node,Tree*t);

	void update_Aug(Node*node,Tree*t);
	kdAugNode_MIN*createNode();
};

class kdLinearAugNode : public kdAugNode
{
public:
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;

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
	kdTree(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat);

	void getNode_Boundary(kdNode*node);
	double obtain_SplitValue(kdNode*node,int split_Dim);
	void KD_Tree_Recur(kdNode*node,int split_Dim);
	void build_kdTree(SVM_stat& stat);
	void initTree_alpha(kdNode*node);

	void updateAugment(kdNode*node);

	//Used for testing the correctness of the tree
	void lookUp_KDTree();
	void lookUp_KDTree_Recur(kdNode*node);

	//added this function for testing on 8/3/2018 
	int select_split_Dim(kdNode*node);
};

class kdTree_adv : public kdTree
{
public:
	//Member functions
	kdTree_adv(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat)
		: kdTree(dim,dataMatrix,alphaArray,leafCapacity,stat){}
	void KD_Tree_adv_Recur(kdNode*node,int split_Dim);
	void build_kdTree_adv(SVM_stat& stat);
	void merge_Boundary(kdNode*node,double**boundary_1,double**boundary_2);
};

#endif