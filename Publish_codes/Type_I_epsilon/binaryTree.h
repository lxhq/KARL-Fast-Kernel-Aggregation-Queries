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

	double LB(double*q,int dim,KDE_stat& stat); 
	double UB(double*q,int dim,KDE_stat& stat); 

	void update_Aug(Node*node,Tree*t); 
	virtual void update_Aug(double**dataMatrix,int dim,KDE_stat& stat); 

	//finish
	virtual void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,int dim,KDE_stat& stat); //GPA

	void initNode(int dim);

	void Aug_Incr(double**dataMatrix,int dim,KDE_stat& stat,int id);
	void Aug_Decr(double**dataMatrix,int dim,KDE_stat& stat,int id);
	void update_Radius(double**dataMatrix,int dim);
	void assign(binaryNode*bNode_delta,int dim);

	binaryNode*createNode();
};

class ballNode : public binaryNode
{
public:
	double**boundary;
	//LB and UB functions declaration
	double UB(double*q,int dim,KDE_stat& stat);

	~ballNode(){}
	void update_Aug(double**dataMatrix,int dim,KDE_stat& stat);
	//add GPA function
	void GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,int dim,KDE_stat& stat);
	void updateBoundary(double**dataMatrix,KDE_stat& stat,int dim);
	ballNode*createNode();
};

class ballNode_SOTA : public ballNode
{
public:
	~ballNode_SOTA(){}
	double LB(double*q,int dim,KDE_stat& stat);
	double UB(double*q,int dim,KDE_stat& stat);
	ballNode_SOTA*createNode();
};

class binaryTree : public Tree
{
public:
	binaryTree(int dim,double**dataMatrix,int leafCapacity,KDE_stat& stat);
	
	void build_BinaryTreeRecur(binaryNode*node);
	void build_BinaryTree();
};

#endif