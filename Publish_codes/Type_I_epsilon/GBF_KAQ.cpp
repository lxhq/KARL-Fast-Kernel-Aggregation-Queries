#include "GBF_KAQ.h"

inline void clearHeap(PQ& pq)
{
	int heapSize=(int)pq.size();
	for(int h=0;h<heapSize;h++)
		pq.pop();
}

void GBF_iter(double*q,Tree& tree,int dim,KDE_stat& stat)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L,U;
	double f_cur;
	double val_R;
	double sq_Euclid;

	Node*rootNode=tree.rootNode;

	L=rootNode->LB(q,dim,stat);
	U=rootNode->UB(q,dim,stat);

	pq_entry.node=rootNode;
	pq_entry.node_L=L;
	pq_entry.node_U=U;
	pq_entry.discrepancy=U-L;

	pq.push(pq_entry);

	while(pq.size()!=0)
	{
		//validation condition
		if(validate_best(L,U,stat.rel_error,val_R)==true)
		{	
			stat.resultValueVector.push_back(val_R);
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
			f_cur=0;
			for(int i=0;i<(int)curNode->idList.size();i++)
			{
				sq_Euclid=0;
				for(int d=0;d<dim;d++)
					sq_Euclid+=(q[d]-tree.dataMatrix[curNode->idList[i]][d])*(q[d]-tree.dataMatrix[curNode->idList[i]][d]);

				f_cur+=exp(-sq_Euclid);
			}

			L=L+f_cur;
			U=U+f_cur;

			continue;
		}

		//Non-Leaf Node
		for(int c=0;c<(int)curNode->childVector.size();c++)
		{
			pq_entry.node_L=curNode->childVector[c]->LB(q,dim,stat);
			pq_entry.node_U=curNode->childVector[c]->UB(q,dim,stat);
			pq_entry.discrepancy=pq_entry.node_U-pq_entry.node_L;
			pq_entry.node=curNode->childVector[c];

			L=L+pq_entry.node_L;
			U=U+pq_entry.node_U;

			pq.push(pq_entry);
		}
	}

	//Case (L=exact=U):
	stat.resultValueVector.push_back(L);
	clearHeap(pq);
}

double computeSqNorm(double*q,int dim)
{
	double sqNorm=0;
	for(int d=0;d<dim;d++)
		sqNorm+=q[d]*q[d];

	return sqNorm;
}

void KAQ_Algorithm(double**queryMatrix,double**dataMatrix,int qNum,int dim,int leafCapacity,int internalCapacity,int method,KDE_stat& stat,char*bulkLoad_TreeName)
{
	#ifndef C_PLUSPLUS11_CLOCK
		clock_t start_s;
		clock_t end_s;
	#endif

	double online_Time;
	kdTree kd_Tree(dim,dataMatrix,leafCapacity,stat);
	mTree m_Tree(dim,dataMatrix,internalCapacity,leafCapacity,stat);
	binaryTree binary_Tree(dim,dataMatrix,leafCapacity,stat);

	//init rootNode
	if(method==1)//tKDC kd-tree + LB_MBR and UB_MBR
		kd_Tree.rootNode=new kdNode();
	if(method==2)//kd-tree + LB_Tri and UB_Tri
		kd_Tree.rootNode=new kdAugNode();
	if(method==3)//kd-tree + Linear Bound
		kd_Tree.rootNode=new kdLinearAugNode();
	if(method==4)//m-tree + triangle inequality
		m_Tree.rootNode=new mNode();
	if(method==5)//m-tree + LB_MBR and UB_MBR
		m_Tree.rootNode=new mAugNode();
	if(method==6)//m-tree + LinearApprox (Tri functions in ell/u)
		m_Tree.rootNode=new mLinearAugNode();
	if(method==7)//m-tree + LinearApprox (MBR functions in ell/u)
		m_Tree.rootNode=new mLinearAug_RectNode();
	if(method==11)//ball-tree + LB_MBR and UB_MBR
		binary_Tree.rootNode=new ballNode_SOTA();
	if(method==12)//ball-tree + LinearApprox
		binary_Tree.rootNode=new ballNode();

	//build tree Preprocessing + create augment tree
	if(method>=1 && method<=3)
	{
		kd_Tree.build_kdTree(stat);
		kd_Tree.updateAugment((kdNode*)kd_Tree.rootNode);
	}

	if(method>=4 && method<=7)
		m_Tree.load_Tree(bulkLoad_TreeName);

	if(method==11 || method==12)
		binary_Tree.build_BinaryTree();

	//no-index algorithm preprocessing
	double**boundary;
	double*center;
	double radius;
	double*a_G;
	double S_G;

	if(method==8)
		pre_Compute_sequential_MBR(dataMatrix,dim,boundary,stat);
	if(method==9)
		pre_Compute_sequential_Delta(dataMatrix,dim,center,radius,stat);
	if(method==10)
		pre_Compute_sequential(dataMatrix,dim,boundary,a_G,S_G,stat);

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		start_s=clock();
	#endif

	//method=0 Linear Scan
	//method=1 kd-tree (LB_MBR, UB_MBR)
	//method=2 kd-tree (LB_Tri, UB_Tri)
	//method=3 kd-tree (LB_(T,t^*), UB_C)
	//method=4 bulk-loading m-tree (LB-tri, UB-tri)
	//method=5 bulk-loading m-tree (LB_MBR, UB_MBR)
	//method=6 bulk-loading m-tree (LB_(T,t^*), UB_C) Tri functions in ell,u
	//method=7 bulk-loading m-tree (LB_(T,t^*), UB_C) MBR functions in ell,u
	//method=8 Sequential Scan with MBR bound functions
	//method=9 Sequential Scan with Delta bound functions
	//method=10 Sequential Scan with our bound functions
	//method=11 ball-tree SOTA
	//method=12 ball-tree KARL
	for(int q=0;q<qNum;q++)
	{
		switch(method)
		{
			case 0:
				SS_iter(queryMatrix[q],dataMatrix,dim,stat);
				break;
			case 1:
			case 2:
				GBF_iter(queryMatrix[q],kd_Tree,dim,stat);
				break;
			case 3:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],kd_Tree,dim,stat);
				break;
			case 4:
			case 5:
				GBF_iter(queryMatrix[q],m_Tree,dim,stat);
				break;
			case 6:
			case 7:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],m_Tree,dim,stat);
				break;
			case 8:
				SS_MBR(queryMatrix[q],dataMatrix,dim,boundary,stat);
				break;
			case 9:
				SS_Delta(queryMatrix[q],dataMatrix,dim,center,radius,stat);
				break;
			case 10:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				SS_linear(queryMatrix[q],dataMatrix,dim,boundary,a_G,S_G,stat);
				break;
			case 11:
				GBF_iter(queryMatrix[q],binary_Tree,dim,stat);
				break;
			case 12:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],binary_Tree,dim,stat);
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

	cout<<"Method "<<method<<": "<<((double)qNum/online_Time)<<" Queries/sec"<<endl;
	//cout<<"pruning ratio: "<<((double)stat.pruneCount)/((double)qNum)<<endl;

}