#include "GBF_KAQ.h"
#include <chrono>
#include <iomanip>

inline void clearHeap(PQ& pq)
{
	int heapSize=(int)pq.size();
	for(int h=0;h<heapSize;h++)
		pq.pop();
}

static double elapsedSeconds(std::chrono::high_resolution_clock::time_point start,
							 std::chrono::high_resolution_clock::time_point end)
{
	return std::chrono::duration<double>(end-start).count();
}

static double percentileLongLong(vector<long long> values,double percentile)
{
	if(values.empty())
		return 0.0;

	sort(values.begin(),values.end());
	int index=(int)ceil(percentile*((double)values.size()))-1;
	if(index<0)
		index=0;
	if(index>=(int)values.size())
		index=(int)values.size()-1;

	return (double)values[index];
}

static void printProfile(int method,KDE_stat& stat)
{
	if(stat.profile_queries==0)
		return;

	double query_count=(double)stat.profile_queries;
	double avg_exact_points=((double)stat.profile_exact_points)/query_count;
	double exact_ratio=stat.n==0 ? 0.0 : avg_exact_points/((double)stat.n);

	cout<<fixed<<setprecision(6);
	cout<<"[PROFILE][M"<<method<<"] queries="<<stat.profile_queries<<", dataset_n="<<stat.n<<endl;
	cout<<"[PROFILE][M"<<method<<"] avg_nodes_processed_per_query="
		<<((double)stat.profile_nodes_processed)/query_count
		<<", p95_nodes_processed_per_query="
		<<percentileLongLong(stat.profile_nodes_per_query,0.95)<<endl;
	cout<<"[PROFILE][M"<<method<<"] avg_leaf_nodes_per_query="
		<<((double)stat.profile_leaf_nodes)/query_count
		<<", p95_leaf_nodes_per_query="
		<<percentileLongLong(stat.profile_leaf_nodes_per_query,0.95)<<endl;
	cout<<"[PROFILE][M"<<method<<"] avg_exact_points_per_query="
		<<avg_exact_points
		<<", p95_exact_points_per_query="
		<<percentileLongLong(stat.profile_exact_points_per_query,0.95)
		<<", exact_point_eval_ratio_vs_full_scan="<<exact_ratio<<endl;
	cout<<"[PROFILE][M"<<method<<"] avg_bound_calls_per_query="
		<<((double)stat.profile_bound_calls)/query_count
		<<", avg_validate_checks_per_query="
		<<((double)stat.profile_validate_checks)/query_count<<endl;
	cout<<"[PROFILE][M"<<method<<"] avg_heap_pushes_per_query="
		<<((double)stat.profile_heap_pushes)/query_count
		<<", avg_heap_pops_per_query="
		<<((double)stat.profile_heap_pops)/query_count<<endl;
	cout<<"[PROFILE][M"<<method<<"] validate_success_queries="
		<<stat.profile_validate_success_queries
		<<", exact_finish_queries="<<stat.profile_exact_finish_queries<<endl;
	cout<<"[PROFILE][M"<<method<<"] time_validate_sec="
		<<stat.profile_time_validate_sec
		<<", time_bound_sec="<<stat.profile_time_bound_sec
		<<", time_leaf_eval_sec="<<stat.profile_time_leaf_eval_sec<<endl;
	if(stat.profile_residual_discards>0 || stat.profile_residual_budget>0)
	{
		cout<<"[PROFILE][M"<<method<<"] avg_residual_discards_per_query="
			<<((double)stat.profile_residual_discards)/query_count
			<<", avg_final_residual_budget_per_query="
			<<stat.profile_residual_budget/query_count<<endl;
	}
}

static double residualBudgetLimit(double lower,double rel_error,double budget_fraction)
{
	if(lower<=0 || rel_error>=1.0)
		return 0.0;

	return budget_fraction*(2.0*rel_error*lower)/(1.0-rel_error);
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
	long long query_nodes_processed=0;
	long long query_leaf_nodes=0;
	long long query_exact_points=0;
	long long query_bound_calls=0;
	long long query_validate_checks=0;
	long long query_heap_pushes=0;
	long long query_heap_pops=0;

	auto recordQueryProfile=[&](bool validateSuccess,bool exactFinish)
	{
		stat.profile_queries++;
		stat.profile_nodes_processed+=query_nodes_processed;
		stat.profile_leaf_nodes+=query_leaf_nodes;
		stat.profile_exact_points+=query_exact_points;
		stat.profile_bound_calls+=query_bound_calls;
		stat.profile_validate_checks+=query_validate_checks;
		stat.profile_heap_pushes+=query_heap_pushes;
		stat.profile_heap_pops+=query_heap_pops;
		stat.profile_nodes_per_query.push_back(query_nodes_processed);
		stat.profile_leaf_nodes_per_query.push_back(query_leaf_nodes);
		stat.profile_exact_points_per_query.push_back(query_exact_points);
		if(validateSuccess)
			stat.profile_validate_success_queries++;
		if(exactFinish)
			stat.profile_exact_finish_queries++;
	};

	Node*rootNode=tree.rootNode;

	auto bound_start=std::chrono::high_resolution_clock::now();
	L=rootNode->LB(q,dim,stat);
	U=rootNode->UB(q,dim,stat);
	auto bound_end=std::chrono::high_resolution_clock::now();
	stat.profile_time_bound_sec+=elapsedSeconds(bound_start,bound_end);
	query_bound_calls+=2;

	pq_entry.node=rootNode;
	pq_entry.node_L=L;
	pq_entry.node_U=U;
	pq_entry.discrepancy=U-L;
	pq.push(pq_entry);
	query_heap_pushes++;

	while(pq.size()!=0)
	{
		auto validate_start=std::chrono::high_resolution_clock::now();
		bool validated=validate_best(L,U,stat.rel_error,val_R);
		auto validate_end=std::chrono::high_resolution_clock::now();
		stat.profile_time_validate_sec+=elapsedSeconds(validate_start,validate_end);
		query_validate_checks++;

		if(validated==true)
		{
			stat.resultValueVector.push_back(val_R);
			recordQueryProfile(true,false);
			clearHeap(pq);
			return;
		}

		pq_entry=pq.top();
		pq.pop();
		query_heap_pops++;
		query_nodes_processed++;

		L=L-pq_entry.node_L;
		U=U-pq_entry.node_U;

		curNode=pq_entry.node;

		//leaf Node
		if((int)curNode->idList.size()<=tree.leafCapacity)
		{
			f_cur=0;
			query_leaf_nodes++;
			query_exact_points+=(long long)curNode->idList.size();
			auto leaf_start=std::chrono::high_resolution_clock::now();
			for(int i=0;i<(int)curNode->idList.size();i++)
			{
				sq_Euclid=0;
				for(int d=0;d<dim;d++)
					sq_Euclid+=(q[d]-tree.dataMatrix[curNode->idList[i]][d])*(q[d]-tree.dataMatrix[curNode->idList[i]][d]);

				f_cur+=exp(-sq_Euclid);
			}
			auto leaf_end=std::chrono::high_resolution_clock::now();
			stat.profile_time_leaf_eval_sec+=elapsedSeconds(leaf_start,leaf_end);

			L=L+f_cur;
			U=U+f_cur;

			continue;
		}

		//Non-Leaf Node
		for(int c=0;c<(int)curNode->childVector.size();c++)
		{
			auto child_bound_start=std::chrono::high_resolution_clock::now();
			pq_entry.node_L=curNode->childVector[c]->LB(q,dim,stat);
			pq_entry.node_U=curNode->childVector[c]->UB(q,dim,stat);
			auto child_bound_end=std::chrono::high_resolution_clock::now();
			stat.profile_time_bound_sec+=elapsedSeconds(child_bound_start,child_bound_end);
			query_bound_calls+=2;
			pq_entry.discrepancy=pq_entry.node_U-pq_entry.node_L;
			pq_entry.node=curNode->childVector[c];

			L=L+pq_entry.node_L;
			U=U+pq_entry.node_U;

			pq.push(pq_entry);
			query_heap_pushes++;
		}
	}

	//Case (L=exact=U):
	stat.resultValueVector.push_back(L);
	recordQueryProfile(false,true);
	clearHeap(pq);
}

void GBF_iter_residual_budget(double*q,Tree& tree,int dim,KDE_stat& stat,double budget_fraction)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L,U;
	double residual_B=0;
	double f_cur;
	double val_R;
	double sq_Euclid;
	vector<pqNode> residual_nodes;
	long long query_nodes_processed=0;
	long long query_leaf_nodes=0;
	long long query_exact_points=0;
	long long query_bound_calls=0;
	long long query_validate_checks=0;
	long long query_heap_pushes=0;
	long long query_heap_pops=0;
	long long query_residual_discards=0;

	auto recordQueryProfile=[&](bool validateSuccess,bool exactFinish)
	{
		stat.profile_queries++;
		stat.profile_nodes_processed+=query_nodes_processed;
		stat.profile_leaf_nodes+=query_leaf_nodes;
		stat.profile_exact_points+=query_exact_points;
		stat.profile_bound_calls+=query_bound_calls;
		stat.profile_validate_checks+=query_validate_checks;
		stat.profile_heap_pushes+=query_heap_pushes;
		stat.profile_heap_pops+=query_heap_pops;
		stat.profile_residual_discards+=query_residual_discards;
		stat.profile_residual_budget+=residual_B;
		stat.profile_nodes_per_query.push_back(query_nodes_processed);
		stat.profile_leaf_nodes_per_query.push_back(query_leaf_nodes);
		stat.profile_exact_points_per_query.push_back(query_exact_points);
		if(validateSuccess)
			stat.profile_validate_success_queries++;
		if(exactFinish)
			stat.profile_exact_finish_queries++;
	};

	Node*rootNode=tree.rootNode;

	auto bound_start=std::chrono::high_resolution_clock::now();
	L=rootNode->LB(q,dim,stat);
	U=rootNode->UB(q,dim,stat);
	auto bound_end=std::chrono::high_resolution_clock::now();
	stat.profile_time_bound_sec+=elapsedSeconds(bound_start,bound_end);
	query_bound_calls+=2;

	pq_entry.node=rootNode;
	pq_entry.node_L=L;
	pq_entry.node_U=U;
	pq_entry.discrepancy=U-L;
	pq.push(pq_entry);
	query_heap_pushes++;

	while(true)
	{
		while(pq.size()!=0)
		{
			auto validate_start=std::chrono::high_resolution_clock::now();
			bool validated=validate_best(L,U+residual_B,stat.rel_error,val_R);
			auto validate_end=std::chrono::high_resolution_clock::now();
			stat.profile_time_validate_sec+=elapsedSeconds(validate_start,validate_end);
			query_validate_checks++;

			if(validated==true)
			{
				stat.resultValueVector.push_back(val_R);
				recordQueryProfile(true,false);
				clearHeap(pq);
				return;
			}

			pq_entry=pq.top();
			pq.pop();
			query_heap_pops++;
			query_nodes_processed++;

			L=L-pq_entry.node_L;
			U=U-pq_entry.node_U;

			curNode=pq_entry.node;

			if((int)curNode->idList.size()<=tree.leafCapacity)
			{
				f_cur=0;
				query_leaf_nodes++;
				query_exact_points+=(long long)curNode->idList.size();
				auto leaf_start=std::chrono::high_resolution_clock::now();
				for(int i=0;i<(int)curNode->idList.size();i++)
				{
					sq_Euclid=0;
					for(int d=0;d<dim;d++)
						sq_Euclid+=(q[d]-tree.dataMatrix[curNode->idList[i]][d])*(q[d]-tree.dataMatrix[curNode->idList[i]][d]);

					f_cur+=exp(-sq_Euclid);
				}
				auto leaf_end=std::chrono::high_resolution_clock::now();
				stat.profile_time_leaf_eval_sec+=elapsedSeconds(leaf_start,leaf_end);

				L=L+f_cur;
				U=U+f_cur;

				continue;
			}

			for(int c=0;c<(int)curNode->childVector.size();c++)
			{
				auto child_bound_start=std::chrono::high_resolution_clock::now();
				pq_entry.node_L=curNode->childVector[c]->LB(q,dim,stat);
				pq_entry.node_U=curNode->childVector[c]->UB(q,dim,stat);
				auto child_bound_end=std::chrono::high_resolution_clock::now();
				stat.profile_time_bound_sec+=elapsedSeconds(child_bound_start,child_bound_end);
				query_bound_calls+=2;

				pq_entry.discrepancy=pq_entry.node_U-pq_entry.node_L;
				pq_entry.node=curNode->childVector[c];

				double budget_limit=residualBudgetLimit(L,stat.rel_error,budget_fraction);
				if(residual_B+pq_entry.node_U<=budget_limit)
				{
					residual_B+=pq_entry.node_U;
					residual_nodes.push_back(pq_entry);
					query_residual_discards++;
					continue;
				}

				L=L+pq_entry.node_L;
				U=U+pq_entry.node_U;

				pq.push(pq_entry);
				query_heap_pushes++;
			}
		}

		if(validate_best(L,U+residual_B,stat.rel_error,val_R)==true)
		{
			stat.resultValueVector.push_back(val_R);
			recordQueryProfile(true,false);
			clearHeap(pq);
			return;
		}

		if(residual_nodes.empty())
			break;

		for(int i=0;i<(int)residual_nodes.size();i++)
		{
			L+=residual_nodes[i].node_L;
			U+=residual_nodes[i].node_U;
			pq.push(residual_nodes[i]);
			query_heap_pushes++;
		}
		residual_nodes.clear();
		residual_B=0;
	}

	stat.resultValueVector.push_back(L);
	recordQueryProfile(false,true);
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
	if(method==13)//kd-tree + Anchor-factorized Linear Bound
		kd_Tree.rootNode=new kdAnchorAugNode();
	if(method==14)//kd-tree + Anchor-factorized Linear Bound only
		kd_Tree.rootNode=new kdAnchorOnlyAugNode();
	if(method==15)//kd-tree + adaptive selective combined KARL/anchor bound
		kd_Tree.rootNode=new kdAdaptiveCombinedAnchorNode();
	if(method==16)//kd-tree + adaptive selective KARL-or-anchor bound
		kd_Tree.rootNode=new kdAdaptiveSelectAnchorNode();
	if(method==17 || method==18)//kd-tree + original KARL bound with residual-budget pruning
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
	if((method>=1 && method<=3) || (method>=13 && method<=18))
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
	//method=13 kd-tree anchor-factorized KARL
	//method=14 kd-tree anchor-factorized KARL without original-bound combination
	//method=15 kd-tree adaptive selective combined KARL/anchor bound
	//method=16 kd-tree adaptive selective KARL-or-anchor bound
	//method=17 kd-tree original KARL with conservative residual-budget pruning
	//method=18 kd-tree original KARL with aggressive residual-budget pruning
	auto start_time = std::chrono::high_resolution_clock::now();
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
			case 13:
			case 14:
			case 15:
			case 16:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],kd_Tree,dim,stat);
				break;
			case 17:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter_residual_budget(queryMatrix[q],kd_Tree,dim,stat,0.5);
				break;
			case 18:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter_residual_budget(queryMatrix[q],kd_Tree,dim,stat,1.0);
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
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
	std::cout << "Total time: " << elapsed.count() << " seconds" << std::endl;
	#ifdef C_PLUSPLUS11_CLOCK
		auto end_s=chrono::high_resolution_clock::now();
		online_Time=(chrono::duration_cast<chrono::nanoseconds>(end_s-start_s).count())/1000000000.0;
	#else
		end_s=clock();
		online_Time=((double)(end_s-start_s))/CLOCKS_PER_SEC;
	#endif

	cout<<"Method "<<method<<": "<<((double)qNum/online_Time)<<" Queries/sec"<<endl;
	printProfile(method,stat);
	//cout<<"pruning ratio: "<<((double)stat.pruneCount)/((double)qNum)<<endl;

}
