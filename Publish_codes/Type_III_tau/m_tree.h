#ifndef M_TREE_H
#define M_TREE_H

#include "Tree.h"

class mNode : public Node
{
public:
	mNode();
	double*O_r_pos;
	double*O_r_neg;
	double radius_pos;
	double radius_neg;
	mNode*parent;

	//facilitates online (sharing) computation
	double temp_obt_dist_pos;
	double temp_obt_dist_neg;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	mNode*createNode(){return new mNode();}
	void update_Aug(Node*node,Tree*t){}

};

class mAugNode : public mNode
{
	double**boundary_pos;
	double**boundary_neg;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);
	void updateBoundary(double**dataMatrix,double*outputArray,SVM_stat& stat,int dim);
	void update_Aug(Node*node,Tree*t);

	mAugNode*createNode(){return new mAugNode();}
};

class mLinearAugNode : public mNode
{
	double*a_G_pos;
	double*a_G_neg;
	double S_G_pos;
	double S_G_neg;

	//facilitates online (sharing) computation
	//double gamma_sum;
	double gamma_sum_pos;
	double gamma_sum_neg;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);
	mLinearAugNode*createNode(){return new mLinearAugNode();}
	void update_Aug(Node*node,Tree*t);
	void update_a_G(double**dataMatrix,double*outputArray,int dim);
	void update_S_G(double**dataMatrix,double*outputArray,int dim);
};

class mTree : public Tree
{
public:
	mTree(int dim,double**dataMatrix,double*outputArray,int internalCapacity,int leafCapacity,SVM_stat& stat);

	double*O_p1;
	double*O_p2;
	double**distMatrix;
	int internalCapacity;

	//bulk-loading
	void build_BL_m_tree_Recur(mNode*node);
	void build_BL_m_tree();
	void update_BL_Augment(mNode*node);
	void sample(mNode*node,vector<int>& sample_idList,int sampleNum);
	void update_info(mNode*node);

	//add write-tree function
	void save_Tree(char*treeFileName);
	void load_Tree(char*treeFileName);
};

#endif