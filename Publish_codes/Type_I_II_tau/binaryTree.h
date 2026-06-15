#ifndef BINARYTREE_H
#define BINARYTREE_H

#include "Tree.h"

class binaryNode : public Node
{
public:
	virtual ~binaryNode(){}

	//center and radius are for the UB
	double*center;
	double radius;

	//These variables are for the LB
	double*a_G;
	double S_G;

	//facilitates online (sharing) computation
	double gamma_sum;
	double temp_obt_dist;

	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);

	//void update_linearAugInfo(binaryNode*node,Tree*t);
	void update_Aug(Node*node,Tree*t);
	virtual void update_Aug(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);

	virtual void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat); //GPA

	void initNode(int dim);

	void Aug_Incr(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat,int id/*,double& update_Radius*/);
	void Aug_Decr(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat,int id/*,double& update_Radius*/);
	void update_Radius(double**dataMatrix,double*alphaArray,int dim);
	void assign(binaryNode*bNode_delta,int dim);

	binaryNode*createNode();
};

class normNode  : public binaryNode
{
	~normNode(){}
	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	normNode*createNode();
};

class farNode : public binaryNode
{
	~farNode(){}
	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	farNode*createNode();
};

class ipNode : public binaryNode
{
	~ipNode(){}
	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	ipNode*createNode();
};

class ballNode : public binaryNode
{
public:
	double**boundary;
	//LB and UB functions declaration
	double UB(double*q,int dim,SVM_stat& stat);

	~ballNode(){}
	void update_Aug(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	//add GPA function
	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	void updateBoundary(double**dataMatrix,SVM_stat& stat,int dim);
	ballNode*createNode();
};

class ballNode_SOTA : public ballNode
{
public:
	~ballNode_SOTA(){}
	double LB(double*q,int dim,SVM_stat& stat);
	double UB(double*q,int dim,SVM_stat& stat);
	ballNode_SOTA*createNode();
};

class oracleNode : public binaryNode
{
public:
	~oracleNode(){}
	double outerRadius;
	double innerRadius;
	//double*ring_center;

	double UB(double*q,int dim,SVM_stat& stat);

	void updateOracleNode(double*ext,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	void updateOracle_LeafModel(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat); //(8) update model //(9) update x_space

	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat);
	oracleNode*createNode();
};

class oracle_iNode : public oracleNode
{
	double UB(double*q,int dim,SVM_stat& stat);
	oracle_iNode*createNode();
};

class binaryTree : public Tree
{
public:
	binaryTree(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat);
	
	void build_BinaryTreeRecur(binaryNode*node);
	void build_BinaryTree();

	//Used for Oracle-node
	void update_OracleRecur(double*ext,oracleNode*node);
	void update_Oracle(double*ext);

	void lookUp_BinaryTree();
	void lookUp_BinaryTree_Recur(binaryNode*node);
};

class ringTree : public binaryTree
{
public:
	ringTree(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat) 
		: binaryTree(dim,dataMatrix,alphaArray,leafCapacity,stat){}
	double*ring_center;
	void build_RingTree();
	void obtain_RingCenter();
};

#endif