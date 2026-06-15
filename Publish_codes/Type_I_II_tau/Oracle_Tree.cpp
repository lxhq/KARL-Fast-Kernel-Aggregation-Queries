#include "Oracle_Tree.h"

double runOracleTree(double*q,svm_node*svm_q,binaryTree& oracle_Tree)
{
	vector< orderEntry > orderList;
	//reorder the dataMatrix
	oracle_Reorder(q,oracle_Tree.dataMatrix,orderList,oracle_Tree.dim,oracle_Tree.stat);
	reordering(oracle_Tree.dataMatrix,oracle_Tree.alphaArray,orderList,oracle_Tree.stat,oracle_Tree.dim);
	//update the binary tree
	oracle_Tree.update_Oracle(q);

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_s=clock();
	#endif

	GBF_iter(q,svm_q,oracle_Tree,oracle_Tree.dim,oracle_Tree.stat);

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_s=chrono::high_resolution_clock::now();
		double online_Time=(chrono::duration_cast<chrono::nanoseconds>(end_s-start_s).count());
	#else
		clock_t end_s=clock();
		double online_Time=((double)(end_s-start_s))/CLOCKS_PER_SEC;
	#endif

	orderList.clear();

	return online_Time;
}

void runOracleTree_mless(double*q,svm_node*svm_q,double**& PP,double*& reorder_alphaArray,binaryTree& oracle_Tree,double& online_Time)
{
	vector< orderEntry > orderList;
	//double online_Time;

	//reorder the dataMatrix
	oracle_Reorder(q,oracle_Tree.dataMatrix,orderList,oracle_Tree.dim,oracle_Tree.stat);
	reordering_mless(oracle_Tree.dataMatrix,PP,oracle_Tree.alphaArray,reorder_alphaArray,orderList,oracle_Tree.stat,oracle_Tree.dim);
	orderList.clear();

	//update the binary tree
	oracle_Tree.update_Oracle(q);

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_s=clock();
	#endif

	GBF_iter(q,svm_q,oracle_Tree,oracle_Tree.dim,oracle_Tree.stat);

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_s=chrono::high_resolution_clock::now();
		online_Time=(chrono::duration_cast<chrono::nanoseconds>(end_s-start_s).count());
	#else
		clock_t end_s=clock();
		online_Time=((double)(end_s-start_s))/CLOCKS_PER_SEC;
	#endif

	//return online_Time;
}