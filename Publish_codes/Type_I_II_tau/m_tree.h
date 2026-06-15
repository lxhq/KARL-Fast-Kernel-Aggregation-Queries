#ifndef M_TREE_H
#define M_TREE_H

#include "Tree.h"

class mNode : public Node
{
public:
	mNode();
	double*O_r;
	double radius;
	mNode*parent;

	//facilitates online (sharing) computation
	double temp_obt_dist;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	mNode*createNode(){return new mNode();}
	void update_Aug(Node*node,Tree*t){}

};

class mAugNode : public mNode
{
public:
	double**boundary;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);
	void updateBoundary(double**dataMatrix,SVM_stat& stat,int dim);
	void update_Aug(Node*node,Tree*t);

	mAugNode*createNode(){return new mAugNode();}
};

//class mLinearAugNode : public mNode
class mLinearAugNode : public mAugNode
{
public:
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);
	mLinearAugNode*createNode(){return new mLinearAugNode();}
	void update_Aug(Node*node,Tree*t);
	void update_a_G(double**dataMatrix,double*alphaArray,int dim);
	void update_S_G(double**dataMatrix,double*alphaArray,int dim);
};

class mTree : public Tree
{
public:
	mTree(int dim,double**dataMatrix,double*alphaArray,int internalCapacity,int leafCapacity,SVM_stat& stat);

	double*O_p1;
	double*O_p2;
	double**distMatrix;
	int internalCapacity;

	//Original
	void insert(mNode*node,int id);
	void split(mNode*node,bool isLeaf);
	void promote(mNode*node,bool isLeaf);
	void partition(mNode*node,mNode*node_dash,bool isLeaf);
	void build_m_tree();
	void lookUp_m_tree();
	void lookUp_m_tree_Recur(mNode*node);
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

#endif