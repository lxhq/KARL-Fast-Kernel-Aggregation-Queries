#include "GBF_KC_Binary.h"

inline void clearHeap(PQ& pq)
{
	int heapSize=(int)pq.size();
	for(int h=0;h<heapSize;h++)
		pq.pop();
}

double computeSqNorm(double*q,int dim)
{
	double sqNorm=0;
	for(int d=0;d<dim;d++)
		sqNorm+=q[d]*q[d];

	return sqNorm;
}

//The tree can be any tree
void GBF_iter(double*q,svm_node*svm_q,Tree& tree,int dim,SVM_stat& stat)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L,U;
	double f_cur;

	Node*rootNode=tree.rootNode;

	#ifdef PRUNE_STATS
		stat.numBound+=2.0;
	#endif

	#ifdef COMPONENT_CLOCK
		#ifdef C_PLUSPLUS11_CLOCK
			auto start_bound_s=chrono::high_resolution_clock::now();
		#endif
	#endif

	L=rootNode->LB(q,dim,stat);
	U=rootNode->UB(q,dim,stat);

	#ifdef COMPONENT_CLOCK
		#ifdef C_PLUSPLUS11_CLOCK
			auto end_bound_s=chrono::high_resolution_clock::now();
			stat.BC_Time+=(chrono::duration_cast<chrono::nanoseconds>(end_bound_s-start_bound_s).count())/1000000000.0;
		#endif
	#endif

	pq_entry.node=rootNode;
	pq_entry.node_L=L;
	pq_entry.node_U=U;
	pq_entry.discrepancy=U-L;

	pq.push(pq_entry);

	while(pq.size()!=0)
	{
		if(L>=stat.rho)
		{
			stat.class_resultVector.push_back(1);
			stat.pruneCount++;
			clearHeap(pq);

			return;
		}
		if(U<stat.rho)
		{
			stat.class_resultVector.push_back(-1);
			stat.pruneCount++;
			clearHeap(pq);

			return;
		}

		pq_entry=pq.top();
		pq.pop();

		L=L-pq_entry.node_L;
		U=U-pq_entry.node_U;

		curNode=pq_entry.node;

		//leaf Node
		if((int)curNode->idList.size()<=tree.leafCapacity)
		{
			#ifdef PRUNE_STATS
				stat.numExact+=(double)curNode->idList.size();
			#endif

			#ifdef COMPONENT_CLOCK
				#ifdef C_PLUSPLUS11_CLOCK
					auto start_refine_s=chrono::high_resolution_clock::now();
				#endif
			#endif

			f_cur=callSVM_Refine(svm_q,curNode->model,0,curNode->idList.size());

			#ifdef COMPONENT_CLOCK
				#ifdef C_PLUSPLUS11_CLOCK
					auto end_refine_s=chrono::high_resolution_clock::now();
					stat.EC_Time+=(chrono::duration_cast<chrono::nanoseconds>(end_refine_s-start_refine_s).count())/1000000000.0;
				#endif
			#endif

			L=L+f_cur;
			U=U+f_cur;

			continue;
		}

		#ifdef COMPONENT_CLOCK
			#ifdef C_PLUSPLUS11_CLOCK
				auto start_bound_s=chrono::high_resolution_clock::now();
			#endif
		#endif
		//Non-Leaf Node
		for(int c=0;c<(int)curNode->childVector.size();c++)
		{
			#ifdef PRUNE_STATS
				stat.numBound+=2.0;
			#endif

			pq_entry.node_L=curNode->childVector[c]->LB(q,dim,stat);
			pq_entry.node_U=curNode->childVector[c]->UB(q,dim,stat);
			pq_entry.discrepancy=pq_entry.node_U-pq_entry.node_L;
			pq_entry.node=curNode->childVector[c];

			L=L+pq_entry.node_L;
			U=U+pq_entry.node_U;

			pq.push(pq_entry);
		}

		#ifdef COMPONENT_CLOCK
			#ifdef C_PLUSPLUS11_CLOCK
				auto end_bound_s=chrono::high_resolution_clock::now();
				stat.BC_Time+=(chrono::duration_cast<chrono::nanoseconds>(end_bound_s-start_bound_s).count())/1000000000.0;
			#endif
		#endif
	}

	if(L>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		clearHeap(pq);
		return;
	}
	if(U<stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		clearHeap(pq);
		return;
	}
}

void KC_Algorithm(double**queryMatrix,double**dataMatrix,double*outputArray,int qNum,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,char*bulkLoad_TreeName)
{
	#ifndef C_PLUSPLUS11_CLOCK
		clock_t start_s;
		clock_t end_s;
	#endif

	double online_Time;

	mTree m_Tree(dim,dataMatrix,outputArray,internalCapacity,leafCapacity,stat);
	kdTree kd_Tree(dim,dataMatrix,outputArray,leafCapacity,stat);
	kdTree_adv kd_Tree_adv(dim,dataMatrix,outputArray,leafCapacity,stat);

	//method=8 or 9 or 10 variables
	//**********************//
	double sum_alpha_pos,sum_alpha_neg;
	double*a_G_pos;
	double*a_G_neg;
	double S_G_pos,S_G_neg;
	double*center_pos;
	double*center_neg;
	double radius_pos,radius_neg;
	//**********************//

	//method=9 variables
	double**boundary_pos;
	double**boundary_neg;

	//Used in LibSVM linear scan method
	svm_node**svm_qMatrix;
	svm_node*x_space;
	svm_model*model;

	//It is not used if LibSVM is not used
	if(method!=0)
		convert_2D_qMatrix_to_SVMFormat(queryMatrix,svm_qMatrix,dim,qNum);

	//preprocessing for the methods which require LibSVM
	if(method==7 || method==8 || method==9 || method==10)
	{
		init_model_inMemory(dataMatrix,dim,stat,model,x_space);
		createModel_inMemory(dataMatrix,outputArray,dim,stat,model,x_space);
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto build_time_s=chrono::high_resolution_clock::now();
	#else
		double build_time_s=clock();
	#endif

	if(method==8)
		pre_Compute_Bplus(dataMatrix,outputArray,dim,sum_alpha_pos,sum_alpha_neg,a_G_pos,a_G_neg,S_G_pos,S_G_neg,center_pos,center_neg,radius_pos,radius_neg,stat);
	if(method==9)
		pre_Compute_MBR(dataMatrix,outputArray,dim,sum_alpha_pos,sum_alpha_neg,boundary_pos,boundary_neg,stat);
	if(method==10)
		pre_Compute_Delta(dataMatrix,outputArray,dim,sum_alpha_pos,sum_alpha_neg,center_pos,center_neg,radius_pos,radius_neg,stat);

	//init rootNode
	if(method==1)//method=1 kd-tree with MBR inequality
		kd_Tree.rootNode=new kdNode();
	if(method==3)
		kd_Tree.rootNode=new kdLinearAugNode();
	if(method==4)//method=4 Bulk-loading m-tree with triangle inequality (mNode)
		m_Tree.rootNode=new mNode();
	if(method==5)//method=5 Bulk-loading m-tree with MBR inequality
		m_Tree.rootNode=new mAugNode();
	if(method==6)//method=6 Bulk-loading m-tree with triangle inequality (mLinearAugNode)
		m_Tree.rootNode=new mLinearAugNode();
	if(method==11)
		kd_Tree_adv.rootNode=new kdLinearAugNode_adv();


	if(method==1 || method==3)
	{
		kd_Tree.build_kdTree(stat);
		kd_Tree.updateAugment((kdNode*)kd_Tree.rootNode);
	}
		
	//build tree Preprocessing + create augment tree
	if(method==4 || method==5 || method==6)
		m_Tree.load_Tree(bulkLoad_TreeName);

	if(method==11)
	{
		kd_Tree_adv.build_kdTree_adv(stat);
		kd_Tree_adv.updateAugment((kdNode*)kd_Tree_adv.rootNode);
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto build_time_e=chrono::high_resolution_clock::now();
		double offline_Time=(chrono::duration_cast<chrono::nanoseconds>(build_time_e-build_time_s).count())/1000000000.0;
	#else
		double build_time_e=clock();
		double offline_Time=((double)(build_time_e-build_time_s))/CLOCKS_PER_SEC;
	#endif

	//cout<<"The preprocessing time is: "<<offline_Time<<" sec"<<endl;

	#ifdef PRUNE_STATS
		stat.numBound=0;
		stat.numExact=0;
		stat.boundTime=0;
		stat.exactTime=0;
	#endif

	#ifdef COMPONENT_CLOCK
		stat.EC_Time=0;
		stat.BC_Time=0;
	#endif

	//method=0 Linear Scan
	//method=1 (LB_MBR, UB_MBR + kd_Tree)
	//method=2 (LB_Tri, UB_Tri + kd_Tree)
	//method=3 (LB_T*,UB_C + kd_Tree)
	//method=4 Bulk-loading m-tree with triangle inequality (mNode)
	//method=5 Bulk-loading m-tree with triangle inequality (mAugNode)
	//method=6 Bulk-loading m-tree with triangle inequality (mLinearAugNode)
	//method=7 LibSVM
	//method=8 LibSVM_{B+}
	//method=9 LibSVM_{MBR}
	//method=10 LibSVM_{Delta}
	//method=11 (LB_T*,UB_C + kd_Tree advance version)

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		start_s=clock();
	#endif

	for(int q=0;q<qNum;q++)
	{
		//Online One-Time processing
		switch (method)
		{
			case 0:
				SS_iter(queryMatrix[q],dataMatrix,outputArray,dim,stat);
				break;
			case 1:
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 3:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 4:
			case 5:
				GBF_iter(queryMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
			case 6:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
			case 7:
				stat.class_resultVector.push_back((int)svm_predict(model,svm_qMatrix[q]));
				break;
			case 8:
				LibSVM_Bplus(queryMatrix[q],svm_qMatrix[q],model,dim,sum_alpha_pos,sum_alpha_neg,a_G_pos,a_G_neg,
					S_G_pos,S_G_neg,center_pos,center_neg,radius_pos,radius_neg,stat);
				break;
			case 9:
				LibSVM_MBR(queryMatrix[q],svm_qMatrix[q],model,dim,sum_alpha_pos,sum_alpha_neg,boundary_pos,boundary_neg,stat);
				break;
			case 10:
				LibSVM_Delta(queryMatrix[q],svm_qMatrix[q],model,dim,sum_alpha_pos,sum_alpha_neg,center_pos,center_neg,radius_pos,radius_neg,stat);
				break;
			case 11:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree_adv,dim,stat);
				break;
		}
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_s=chrono::high_resolution_clock::now();
		online_Time=(chrono::duration_cast<chrono::nanoseconds>(end_s-start_s).count())/1000000000.0;
	#else
		end_s=clock();
		online_Time=((double)(end_s-start_s))/CLOCKS_PER_SEC;
	#endif

	cout<<((double)qNum/(offline_Time+online_Time))<<endl;
	//cout<<"Method "<<method<<": "<<((double)qNum/online_Time)<<" Queries/sec"<<endl;
	//cout<<"Method "<<method<<": "<<((double)qNum/(offline_Time+online_Time))<<" Queries/sec (End-to-end)"<<endl;
	//cout<<"pruning ratio: "<<((double)stat.pruneCount)/((double)qNum)<<endl;

	/*#ifdef PRUNE_STATS
		cout<<"numBound: "<<stat.numBound<<endl;
		cout<<"numExact: "<<stat.numExact<<endl;
	#endif

	#ifdef COMPONENT_CLOCK
		cout<<"BC Time (on average): "<<stat.BC_Time/(double)stat.numBound<<endl;
		cout<<"EC Time (on average): "<<stat.EC_Time/(double)stat.numExact<<endl;
	#endif*/
}
