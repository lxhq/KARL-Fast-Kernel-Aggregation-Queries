#include "cache_order.h"

//offline
void init_pivotSet(char*pivotFileName,pivotSet& pSet,int dim)
{
	fstream pivotFile;
	pivotFile.open(pivotFileName);

	if(pivotFile.is_open()==false)
	{
		cout<<"Cannot Open Pivot File!"<<endl;
		exit(1);
	}

	pivotFile>>pSet.num_of_pivots;

	pSet.pivotMatrix=new double*[pSet.num_of_pivots];
	for(int p=0;p<pSet.num_of_pivots;p++)
	{
		pSet.pivotMatrix[p]=new double[dim];

		for(int d=0;d<dim;d++)
			pivotFile>>pSet.pivotMatrix[p][d];
	}

	pivotFile.close();
}

void build_multiple_Tree(vector<binaryTree>& multi_b_Tree,pivotSet& pSet,int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat,int method)
{
	double**PP;
	double*reorder_alphaArray;
	vector<orderEntry> orderList;

	init_reorder_mless(PP,reorder_alphaArray,stat,dim);

	binaryTree binary_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);

	for(int t=0;t<pSet.num_of_pivots;t++)
	{
		multi_b_Tree.push_back(binary_Tree);

		if(method==20)
			multi_b_Tree[t].rootNode=new oracleNode();
		if(method==22)
			multi_b_Tree[t].rootNode=new oracle_iNode();

		oracle_Reorder(pSet.pivotMatrix[t],dataMatrix,orderList,dim,stat);
		reordering_mless(dataMatrix,PP,alphaArray,reorder_alphaArray,orderList,stat,dim);
		orderList.clear();
		
		//build_Tree
		multi_b_Tree[t].build_BinaryTree();
		multi_b_Tree[t].update_Oracle(pSet.pivotMatrix[t]);
	}

	delete_reorder_mless(PP,reorder_alphaArray,stat);
}

//online
int NN(double*q,pivotSet& pSet,int dim)
{
	int tree_num;
	double sq_euclid_dist=0;
	double best_dist=inf;

	for(int p=0;p<pSet.num_of_pivots;p++)
	{
		sq_euclid_dist=0;
		for(int d=0;d<dim;d++)
			sq_euclid_dist+=(q[d]-pSet.pivotMatrix[p][d])*(q[d]-pSet.pivotMatrix[p][d]);

		if(sq_euclid_dist<best_dist)
		{
			best_dist=sq_euclid_dist;
			tree_num=p;
		}
	}

	return tree_num;
}