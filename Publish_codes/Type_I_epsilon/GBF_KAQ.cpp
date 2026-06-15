#include "GBF_KAQ.h"

#ifdef KARL_PROFILE
#include <chrono>
#include <iomanip>
#endif

inline void clearHeap(PQ& pq)
{
	int heapSize=(int)pq.size();
	for(int h=0;h<heapSize;h++)
		pq.pop();
}

#ifdef KARL_PROFILE
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

static const char* stage0Env(const char* name,const char* fallback)
{
	const char* value=getenv(name);
	if(value==NULL || value[0]=='\0')
		return fallback;
	return value;
}

static const char* stage0TreeType(int method)
{
	if(method>=1 && method<=3)
		return "kd_tree";
	if(method==11 || method==12)
		return "ball_tree";
	return "other";
}

static double finalRelativeGap(double L,double U)
{
	double denom=U+L;
	if(fabs(denom)<epsilon)
		return 0.0;
	return (U-L)/denom;
}

static void resetProfile(KDE_stat& stat)
{
	stat.profile_queries=0;
	stat.profile_nodes_processed=0;
	stat.profile_leaf_nodes=0;
	stat.profile_exact_points=0;
	stat.profile_bound_calls=0;
	stat.profile_validate_checks=0;
	stat.profile_heap_pushes=0;
	stat.profile_heap_pops=0;
	stat.profile_validate_success_queries=0;
	stat.profile_exact_finish_queries=0;
	stat.profile_time_validate_sec=0;
	stat.profile_time_bound_sec=0;
	stat.profile_time_leaf_eval_sec=0;
	stat.profile_nodes_per_query.clear();
	stat.profile_leaf_nodes_per_query.clear();
	stat.profile_exact_points_per_query.clear();
	stat.profile_expanded_internal_nodes_per_query.clear();
	stat.profile_accepted_far_nodes_per_query.clear();
	stat.profile_visited_nodes_per_query.clear();
	stat.profile_pq_max_size_per_query.clear();
	stat.profile_certificate_pass_per_query.clear();
	stat.profile_final_relative_gap_per_query.clear();
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
}

static void writeStage0PerQueryCsv(int method,int leafCapacity,int qNum,KDE_stat& stat)
{
	const char* outputDir=getenv("KARL_STAGE0_DIAG_DIR");
	if(outputDir==NULL || outputDir[0]=='\0')
		return;

	const string dataset=stage0Env("KARL_STAGE0_DATASET","unknown");
	const string impl=stage0Env("KARL_STAGE0_IMPL","cpu_karl");
	const string bandwidth=stage0Env("KARL_STAGE0_BANDWIDTH","unknown");
	const string treeType=stage0TreeType(method);
	string outputPath=string(outputDir)+"/"+dataset+".csv";
	ofstream out(outputPath.c_str());
	if(!out)
	{
		cerr<<"[STAGE0] Failed to open per-query CSV: "<<outputPath<<endl;
		return;
	}

	out<<"dataset,impl,tree_type,method,leaf_capacity,bandwidth,query_count,query_id,"
		<<"visited_nodes,expanded_internal_nodes,accepted_far_nodes,visited_leaf_nodes,"
		<<"exact_point_evals,gap_formula,epsilon,certificate_pass,final_relative_gap,"
		<<"pq_max_size,pq_overflow_count\n";

	const string gapFormula="relative=(U-L)/(U+L);pass=U-L<=eps*(U+L);fallback=fabs(L)<1e-9&&fabs(U)<1e-9";
	const int rows=(int)stat.profile_visited_nodes_per_query.size();
	for(int q=0;q<rows;q++)
	{
		out<<dataset<<","<<impl<<","<<treeType<<","<<method<<","<<leafCapacity<<","
			<<bandwidth<<","<<qNum<<","<<q<<","
			<<stat.profile_visited_nodes_per_query[q]<<","
			<<stat.profile_expanded_internal_nodes_per_query[q]<<","
			<<stat.profile_accepted_far_nodes_per_query[q]<<","
			<<stat.profile_leaf_nodes_per_query[q]<<","
			<<stat.profile_exact_points_per_query[q]<<","
			<<"\""<<gapFormula<<"\","
			<<stat.rel_error<<","
			<<stat.profile_certificate_pass_per_query[q]<<","
			<<setprecision(17)<<stat.profile_final_relative_gap_per_query[q]<<","
			<<stat.profile_pq_max_size_per_query[q]<<",0\n";
	}

	cout<<"[STAGE0] Wrote CPU per-query counters: "<<outputPath<<endl;
}
#endif

void GBF_iter(double*q,Tree& tree,int dim,KDE_stat& stat)
{
	static PQ pq;
	pqNode pq_entry;
	Node*curNode;
	double L,U;
	double f_cur;
	double val_R;
	double sq_Euclid;
#ifdef KARL_PROFILE
	long long query_nodes_processed=0;
	long long query_leaf_nodes=0;
	long long query_exact_points=0;
	long long query_expanded_internal_nodes=0;
	long long query_accepted_far_nodes=0;
	long long query_pq_max_size=0;
	long long query_bound_calls=0;
	long long query_validate_checks=0;
	long long query_heap_pushes=0;
	long long query_heap_pops=0;

	auto recordQueryProfile=[&](bool validateSuccess,bool exactFinish,double final_gap)
	{
		long long query_visited_nodes=query_expanded_internal_nodes+query_accepted_far_nodes+query_leaf_nodes;
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
		stat.profile_expanded_internal_nodes_per_query.push_back(query_expanded_internal_nodes);
		stat.profile_accepted_far_nodes_per_query.push_back(query_accepted_far_nodes);
		stat.profile_visited_nodes_per_query.push_back(query_visited_nodes);
		stat.profile_pq_max_size_per_query.push_back(query_pq_max_size);
		stat.profile_certificate_pass_per_query.push_back((validateSuccess || exactFinish) ? 1 : 0);
		stat.profile_final_relative_gap_per_query.push_back(final_gap);
		if(validateSuccess)
			stat.profile_validate_success_queries++;
		if(exactFinish)
			stat.profile_exact_finish_queries++;
	};
#endif

	Node*rootNode=tree.rootNode;

#ifdef KARL_PROFILE
	auto bound_start=std::chrono::high_resolution_clock::now();
#endif
	L=rootNode->LB(q,dim,stat);
	U=rootNode->UB(q,dim,stat);
#ifdef KARL_PROFILE
	auto bound_end=std::chrono::high_resolution_clock::now();
	stat.profile_time_bound_sec+=elapsedSeconds(bound_start,bound_end);
	query_bound_calls+=2;
#endif

	pq_entry.node=rootNode;
	pq_entry.node_L=L;
	pq_entry.node_U=U;
	pq_entry.discrepancy=U-L;

	pq.push(pq_entry);
#ifdef KARL_PROFILE
	query_heap_pushes++;
	query_pq_max_size=max(query_pq_max_size,(long long)pq.size());
#endif

	while(pq.size()!=0)
	{
		//validation condition
#ifdef KARL_PROFILE
		auto validate_start=std::chrono::high_resolution_clock::now();
		bool validated=validate_best(L,U,stat.rel_error,val_R);
		auto validate_end=std::chrono::high_resolution_clock::now();
		stat.profile_time_validate_sec+=elapsedSeconds(validate_start,validate_end);
		query_validate_checks++;
		if(validated==true)
#else
		if(validate_best(L,U,stat.rel_error,val_R)==true)
#endif
		{
			stat.resultValueVector.push_back(val_R);
#ifdef KARL_PROFILE
			query_accepted_far_nodes=(long long)pq.size();
			recordQueryProfile(true,false,finalRelativeGap(L,U));
#endif
			clearHeap(pq);
			return;
		}

		pq_entry=pq.top();
		pq.pop();
#ifdef KARL_PROFILE
		query_heap_pops++;
		query_nodes_processed++;
#endif

		L=L-pq_entry.node_L;
		U=U-pq_entry.node_U;

		curNode=pq_entry.node;

		//leaf Node
		if((int)curNode->idList.size()<=tree.leafCapacity)
		{
			f_cur=0;
#ifdef KARL_PROFILE
			query_leaf_nodes++;
			query_exact_points+=(long long)curNode->idList.size();
			auto leaf_start=std::chrono::high_resolution_clock::now();
#endif
			for(int i=0;i<(int)curNode->idList.size();i++)
			{
				sq_Euclid=0;
				for(int d=0;d<dim;d++)
					sq_Euclid+=(q[d]-tree.dataMatrix[curNode->idList[i]][d])*(q[d]-tree.dataMatrix[curNode->idList[i]][d]);

				f_cur+=exp(-sq_Euclid);
			}
#ifdef KARL_PROFILE
			auto leaf_end=std::chrono::high_resolution_clock::now();
			stat.profile_time_leaf_eval_sec+=elapsedSeconds(leaf_start,leaf_end);
#endif

			L=L+f_cur;
			U=U+f_cur;

			continue;
		}

		//Non-Leaf Node
#ifdef KARL_PROFILE
		query_expanded_internal_nodes++;
#endif
		for(int c=0;c<(int)curNode->childVector.size();c++)
		{
#ifdef KARL_PROFILE
			auto child_bound_start=std::chrono::high_resolution_clock::now();
#endif
			pq_entry.node_L=curNode->childVector[c]->LB(q,dim,stat);
			pq_entry.node_U=curNode->childVector[c]->UB(q,dim,stat);
#ifdef KARL_PROFILE
			auto child_bound_end=std::chrono::high_resolution_clock::now();
			stat.profile_time_bound_sec+=elapsedSeconds(child_bound_start,child_bound_end);
			query_bound_calls+=2;
#endif
			pq_entry.discrepancy=pq_entry.node_U-pq_entry.node_L;
			pq_entry.node=curNode->childVector[c];

			L=L+pq_entry.node_L;
			U=U+pq_entry.node_U;

			pq.push(pq_entry);
#ifdef KARL_PROFILE
			query_heap_pushes++;
			query_pq_max_size=max(query_pq_max_size,(long long)pq.size());
#endif
		}
	}

	//Case (L=exact=U):
	stat.resultValueVector.push_back(L);
#ifdef KARL_PROFILE
	recordQueryProfile(false,true,0.0);
#endif
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

#ifdef KARL_PROFILE
	resetProfile(stat);
#endif

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

#ifdef KARL_PROFILE
	cout<<"Total time: "<<online_Time<<" seconds"<<endl;
#endif
	cout<<"Method "<<method<<": "<<((double)qNum/online_Time)<<" Queries/sec"<<endl;
#ifdef KARL_PROFILE
	printProfile(method,stat);
	writeStage0PerQueryCsv(method,leafCapacity,qNum,stat);
#endif
	//cout<<"pruning ratio: "<<((double)stat.pruneCount)/((double)qNum)<<endl;

}
