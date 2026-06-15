#ifndef M_TREE_H
#define M_TREE_H

#include "Tree.h"

class mNode : public Node
{
public:
	//mNode();
	double*O_r;
	double radius;
	mNode*parent;

	//facilitates online (sharing) computation
	double temp_obt_dist;

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);

	mNode*createNode(){return new mNode();}
	void update_Aug(Node*node,Tree*t){}

};

class mAugNode : public mNode
{
public:
	double**boundary;

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);
	//void updateBoundary(double**dataMatrix,KDE_stat& stat,int dim);
	void update_Aug(Node*node,Tree*t);

	mAugNode*createNode(){return new mAugNode();}
};

class mLinearAugNode : public mNode
{
public:
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;

	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);
	mLinearAugNode*createNode(){return new mLinearAugNode();}
	void update_Aug(Node*node,Tree*t);
	void update_a_G(double**dataMatrix,int dim);
	void update_S_G(double**dataMatrix,int dim);
};

class mLinearAug_RectNode : public mLinearAugNode
{
public:
	double**boundary;
	//Reuse LB function defined in this node (mLinearAugNode)
	double UB(double*q,int dim,KDE_stat& stat);
	void update_Aug(Node*node,Tree*t);
	mLinearAug_RectNode*createNode(){return new mLinearAug_RectNode();}
};

class mTree : public Tree
{
public:
	mTree(int dim,double**dataMatrix,int internalCapacity,int leafCapacity,KDE_stat& stat);

	double*O_p1;
	double*O_p2;
	int internalCapacity;

	void update_Aug(mNode*node); //used in mNode
	void updateAugment(mNode*node); //used in mAugNode and mLinearAugNode
	void update_rootInfo();

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

void updateBoundary(double**dataMatrix,vector<int>& idList,double**& boundary,KDE_stat& stat,int dim);

#endif