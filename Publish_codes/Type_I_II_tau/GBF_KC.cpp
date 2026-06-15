#include "GBF_KC.h"

inline void clearHeap(PQ& pq)
{
	int heapSize=(int)pq.size();
	for(int h=0;h<heapSize;h++)
		pq.pop();
}

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


double computeSqNorm(double*q,int dim)
{
	double sqNorm=0;
	for(int d=0;d<dim;d++)
		sqNorm+=q[d]*q[d];

	return sqNorm;
}

void KC_Algorithm(double**queryMatrix,double**dataMatrix,double*alphaArray,int qNum,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,svm_model*& model,int chunkSize,char*& oracle_TimeFileName,char*pivotFileName,char*bulkLoad_TreeName)
{
	#ifndef C_PLUSPLUS11_CLOCK
		clock_t start_s;
		clock_t end_s;
	#endif

	double online_Time=0;

	kdTree kd_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);
	kdTree_adv kd_Tree_adv(dim,dataMatrix,alphaArray,leafCapacity,stat);
	binaryTree binary_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);
	mTree m_Tree(dim,dataMatrix,alphaArray,internalCapacity,leafCapacity,stat);
	ringTree ring_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);

	//method=9 variables (some of them are used in methods 30,31,32 and 33)
	//**********************//
	double sum_alpha;
	double*a_G;
	double S_G;
	double*center;
	double radius;
	//**********************//

	//method=30/31
	double**boundary;

	//method=11 variables
	//**********************//
	vector<inverted_entry>* inv_Index;
	vector<double> sq_EuclidVector;
	double*p_sqNormArray;
	double**lu_Table;
	//**********************//

	//method=13 variables
	//**********************//
	group_Info*groupArray;
	//**********************//

	//method=100,101,102 variable
	//**********************//
	double bestTime;
	double total_OracleTime=0;
	double**PP;
	double*reorder_alphaArray;
	if(method==102)
		init_reorder_mless(PP,reorder_alphaArray,stat,dim);

	fstream oracle_TimeFile;
	if(method==100 || method==101 || method==102)
	{
		oracle_TimeFile.open(oracle_TimeFileName,ios::out|ios::app);
		if(oracle_TimeFile.is_open()==false)
		{
			cout<<"Cannot Open oracle_TimeFile!"<<endl;
			exit(1);
		}
	}
	//**********************//

	//Cache-based method (preprocessing)
	//**********************//
	int tree_num;
	pivotSet pSet;
	vector<binaryTree> multi_b_Tree;
	if(method==20 || method==22)
	{
		init_pivotSet(pivotFileName,pSet,dim);
		build_multiple_Tree(multi_b_Tree,pSet,dim,dataMatrix,alphaArray,leafCapacity,stat,method);
	}
	//**********************//

	//Used in Reorder Algorithm
	vector<orderEntry> orderList;

	//Used in LibSVM linear scan method
	svm_node**svm_qMatrix;
	convert_2D_qMatrix_to_SVMFormat(queryMatrix,svm_qMatrix,dim,qNum);
	svm_node*x_space;

	double build_Time;
	#ifdef C_PLUSPLUS11_CLOCK
		auto start_b_index_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_b_index_s=clock();
	#endif

	//preprocessing for the methods which require LibSVM
	if(method==9 || method==10 || (method>=12 && method<=19) || method==100 || method==101 || (method>=30 && method<=33))
	{
		init_model_inMemory(dataMatrix,alphaArray,dim,stat,model,x_space);
		if(method==9 || method==10 || method==12 || method==13 || method==16 || (method>=30 && method<=33))
			createModel_inMemory(dataMatrix,alphaArray,dim,stat,model,x_space);
	}

	//testing for method=0
	/*double*dataMatrix_1d;
	if(method==0)
		convert_to_1d(dataMatrix,dataMatrix_1d,dim,stat);*/
	//init rootNode
	if(method==1)//tKDC kd-tree + LB_MBR and UB_MBR
		kd_Tree.rootNode=new kdNode();
	if(method==2)//kd-tree + LB_MBR and UB_MIN (ell_MIN)
		kd_Tree.rootNode=new kdAugNode_MIN();
	if(method==3)//kd-tree + triangle inequality
		kd_Tree.rootNode=new kdAugNode();
	if(method==4)//kd-tree + LinearApprox
		kd_Tree.rootNode=new kdLinearAugNode();
	if(method==5)//binary-tree + LinearApprox + GPA
		binary_Tree.rootNode=new binaryNode();
	if(method==6 || method==27)//m-tree + triangle inequality
		m_Tree.rootNode=new mNode();
	if(method==7 || method==28)//m-tree + LB_MBR and UB_MBR
		m_Tree.rootNode=new mAugNode();
	if(method==8 || method==29)//m-tree + LinearApprox
		m_Tree.rootNode=new mLinearAugNode();
	if(method==23)//ring-tree + LinearApprox
		ring_Tree.rootNode=new oracleNode();
	if(method==24) //binary-tree + LinearApprox + GPA (norm partition)
		binary_Tree.rootNode=new normNode();
	if(method==25)
		binary_Tree.rootNode=new farNode();
	if(method==26)
		binary_Tree.rootNode=new ipNode();
	if(method==34)
		binary_Tree.rootNode=new ballNode();
	if(method==35)
		binary_Tree.rootNode=new ballNode_SOTA();
	if(method==36)
		kd_Tree_adv.rootNode=new kdLinearAugNode_adv();
	if(method==102)//oracle-tree (binary tree)
		binary_Tree.rootNode=new oracleNode();

	//build tree Preprocessing + create augment tree
	if((method>=1 && method<=4) || method==36)
	{
		if(method>=1 && method<=4)
		{
			kd_Tree.build_kdTree(stat);
			kd_Tree.updateAugment((kdNode*)kd_Tree.rootNode);
		}
		if(method==36)
		{
			kd_Tree_adv.build_kdTree_adv(stat);
			kd_Tree_adv.updateAugment((kdLinearAugNode_adv*)kd_Tree_adv.rootNode);
		}
	}
	if(method==5 || method==24 || method==25 || method==26 || method==34 || method==35) //binary tree
		binary_Tree.build_BinaryTree();
		
	if(method>=6 && method<=8)
	{
		m_Tree.build_m_tree();
		m_Tree.updateAugment((mNode*)m_Tree.rootNode);
	}

	if(method==9 || method==10 || method==11)
		pre_Compute_sequential(dataMatrix,alphaArray,dim,sum_alpha,a_G,S_G,center,radius,stat);
	if(method==30 || method==31)
		pre_Compute_sequential_MBR(dataMatrix,alphaArray,dim,sum_alpha,boundary,stat);
	if(method==32 || method==33)
		pre_Compute_sequential_Delta(dataMatrix,alphaArray,dim,sum_alpha,center,radius,stat);

	if(method==11)
	{
		lookUpTable_Creation(lu_Table);
		build_Inverted_Index(dataMatrix,inv_Index,sq_EuclidVector,dim,stat,p_sqNormArray);
	}
	if(method==13 || method==16)
		create_groupArray(dataMatrix,alphaArray,groupArray,chunkSize,dim,stat,true);

	if(method==14 || method==17 || method==19)
	{
		//reorder the dataMatrix and alphaArray
		concentric_Reorder(dataMatrix,orderList,alphaArray,dim,stat);
		reordering(dataMatrix,alphaArray,orderList,stat,dim);
		createModel_inMemory(dataMatrix,alphaArray,dim,stat,model,x_space);
		create_groupArray(dataMatrix,alphaArray,groupArray,chunkSize,dim,stat,true);
	}
	if(method==15 || method==18)
	{
		norm_Reorder(dataMatrix,orderList,dim,stat);
		reordering(dataMatrix,alphaArray,orderList,stat,dim);
		createModel_inMemory(dataMatrix,alphaArray,dim,stat,model,x_space);
		create_groupArray(dataMatrix,alphaArray,groupArray,chunkSize,dim,stat,true);
	}

	if(method==23)
		ring_Tree.build_RingTree();

	if(method==27 || method==28 || method==29)
		m_Tree.load_Tree(bulkLoad_TreeName);

	if(method==102)
		binary_Tree.build_BinaryTree();

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_b_index_s=chrono::high_resolution_clock::now();
		build_Time=(chrono::duration_cast<chrono::nanoseconds>(end_b_index_s-start_b_index_s).count())/1000000000.0;
	#else
		clock_t end_b_index_s=clock();
		build_Time=((double)(end_b_index_s-start_b_index_s))/CLOCKS_PER_SEC;
	#endif

	//cout<<"Preprocessing Time is: "<<build_Time<<" sec"<<endl;

	#ifdef PRUNE_STATS
		stat.numBound=0;
		stat.numExact=0;
		stat.boundTime=0;
		stat.exactTime=0;
		binary_Tree.stat.numBound=0;
		binary_Tree.stat.numExact=0;
		binary_Tree.stat.boundTime=0;
		binary_Tree.stat.exactTime=0;
	#endif

	#ifdef COMPONENT_CLOCK
		stat.EC_Time=0;
		stat.BC_Time=0;
	#endif

	//method=0 Linear Scan
	//method=1 (LB_MBR, UB_MBR + kd_Tree) 
	//method=2 (LB_MBR, UB_MIN + kd_Tree)
	//method=3 (LB_Tri, UB_Tri + kd_Tree)
	//method=4 (LB_T*,UB_C + kd_Tree)
	//method=5 (LB_T*,UB_C + binary_Tree and GPA)
	//method=6 (LB_Tri, UB_Tri + m_Tree)
	//method=7 (LB_MBR, UB_MBR + m_Tree)
	//method=8 (LB_T*,UB_C + m_Tree)
	//method=9 sequential_bound + LB_T*, UB_C (no-tree) via SVM model
	//method=11 sequential_bound + LB_T*, UB_C (no-tree) sparse_eucild via inverted index
	//method=12 LibSVM
	//method=13 incremental classification (IC) via LibSVM
	//method=14 incremental classification (IC) with Concentric Reordering via LibSVM
	//method=15 incremental classification (IC) with Norm Reordering via LibSVM
	//method=16 incremental classification (IC) via LibSVM (multi-step)
	//method=17 incremental classification (IC) with Concentric Reordering via LibSVM (multi-step)
	//method=18 incremental classification (IC) with Norm Reordering via LibSVM (multi-step)
	//method=19 incremental classification (IC) with Concentric Reordering via LibSVM (multi-step) and Statistical Estimation of the chunkSize
	//method=20 cache-based binary-tree (oracle node) with Concentric Reordering
	//method=22 cache-based binary-tree (oracle_i node) with Concentric Reordering
	//method=23 ring-tree with center the centroid point of the dataset
	//method=24 (LB_T*,UB_C + binary_Tree and GPA with norm partition)
	//method=25 (LB_T*,UB_C + binary_Tree and GPA with far partition)
	//method=26 (LB_T*,UB_C + binary_Tree and GPA with ip partition)
	//method=27 Bulk-loading m-tree with triangle inequality (mNode)
	//method=29 Bulk-loading m-tree with triangle inequality (mLinearAugNode)
	//method=30 sequential_bound + LB_{MBR}, UB_{MBR} (no-tree) via SVM model
	//method=31 sequential_bound + LB_{MBR}, UB_{MBR} (no-tree) via SS
	//method=32 sequential_bound + LB_{Delta}, UB_{Delta} (no-tree) via SVM model
	//method=33 sequential_bound + LB_{Delta}, UB_{Delta} (no-tree) via SS
	//method=34 ball-Tree with KARL + MBR-based lb and ub function for UB
	//method=35 ball-Tree with SOTA + MBR-based lb and ub function for UB
	//method=36 kd-Tree with KARL + MBR-based lb and ub function for UB
	//method=37 Weak auto-tuning (choose method=9 or method=36)
	//method=100 Oracle (determine the best number of chunkSize for each query)
	//method=101 Weak Oracle (determine 1/5*total_sv as the chunkSize for all queries, i.e only 5 partitions are used)
	//method=102 Oracle Tree (binary tree)
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
				//SS_iter_1d(queryMatrix[q],dataMatrix_1d,alphaArray,dim,stat);
				SS_iter(queryMatrix[q],dataMatrix,alphaArray,dim,stat);
				break;
			case 1:
			case 2:
			case 3:
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 4:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 5:
			case 24:
			case 25:
			case 26:
			case 34:
			case 35:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],binary_Tree,dim,stat);
				break;
			case 6:
			case 7:
			case 27:
			case 28:
				GBF_iter(queryMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
			case 8:
			case 29:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
			case 9:
			case 10:
				sequential_bound(queryMatrix[q],dataMatrix,alphaArray,svm_qMatrix[q],model,dim,sum_alpha,a_G,S_G,center,radius,stat,method);
				break;
			case 11:
				sequential_bound_sparse(queryMatrix[q],alphaArray,dim,sum_alpha,a_G,S_G,center,radius,stat,inv_Index,sq_EuclidVector,p_sqNormArray,lu_Table);
				break;
			case 12:
				stat.class_resultVector.push_back((int)svm_predict(model,svm_qMatrix[q]));
				break;
			case 13:
			case 14:
			case 15:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				incremental_Classification(queryMatrix[q],dataMatrix,groupArray,alphaArray,chunkSize,dim,stat,model,svm_qMatrix[q]);
				break;
			case 16:
			case 17:
			case 18:
			case 19:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				incremental_Classification_multiStep(queryMatrix[q],dataMatrix,groupArray,alphaArray,chunkSize,dim,stat,model,svm_qMatrix[q]);
				break;
			case 20:
			case 22:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				tree_num=NN(queryMatrix[q],pSet,dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],multi_b_Tree[tree_num],dim,stat);
				break;
			case 23:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],ring_Tree,dim,stat);
				break;
			case 30:
			case 31:
				sequential_bound_MBR(queryMatrix[q],dataMatrix,alphaArray,svm_qMatrix[q],model,dim,sum_alpha,boundary,stat,method);
				break;
			case 32:
			case 33:
				sequential_bound_Delta(queryMatrix[q],dataMatrix,alphaArray,svm_qMatrix[q],model,dim,sum_alpha,center,radius,stat,method);
				break;
			case 36:
			{
				//if(q==147 || q==148 || q==149)
				//	SS_iter(queryMatrix[q],dataMatrix,alphaArray,dim,stat);

				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				GBF_iter(queryMatrix[q],svm_qMatrix[q],kd_Tree_adv,dim,stat);
				break;
			}
			/*case 37: //weak auto-tuning
			{
				int sNum=10;
				method=weak_auto(queryMatrix,dataMatrix,alphaArray,svm_qMatrix,model,sNum,dim,sum_alpha,a_G,S_G,center,radius,stat,leafCapacity,kd_Tree_adv);
				q=2*sNum;
				break;
			}*/
			case 100:
			case 101:
				stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				bestTime=runOracle(queryMatrix[q],dataMatrix,svm_qMatrix[q],model,alphaArray,dim,stat,x_space,method);
				oracle_TimeFile<<bestTime<<endl;
				break;
			case 102:
				binary_Tree.stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
				//bestTime=runOracleTree(queryMatrix[q],svm_qMatrix[q],binary_Tree);
				runOracleTree_mless(queryMatrix[q],svm_qMatrix[q],PP,reorder_alphaArray,binary_Tree,bestTime);
				total_OracleTime=total_OracleTime+bestTime;
				stat.class_resultVector.push_back(binary_Tree.stat.class_resultVector[q]);

				#ifdef PRUNE_STATS
					stat.numBound=binary_Tree.stat.numBound;
					stat.numExact=binary_Tree.stat.numExact;
					stat.pruneCount=binary_Tree.stat.pruneCount;
					stat.boundTime=binary_Tree.stat.boundTime;
					stat.exactTime=binary_Tree.stat.exactTime;
				#endif

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

	if(method!=100 && method!=101 && method!=102)
	{
		cout<</*"Method "<<method<<": "<<*/((double)qNum/online_Time)<<" Queries/sec (Online)"<<endl;
		//cout<<"Method "<<method<<": "<<((double)qNum/(build_Time+online_Time))<<" Queries/sec (End-to-end)"<<endl;
		//cout<<"pruning ratio: "<<((double)stat.pruneCount)/((double)qNum)<<endl;
	}

	if(method==100 || method==101 || method==102)
	{
		cout<<"Method "<<method<<": "<<((double)qNum/(total_OracleTime/1000000000.0))<<" Queries/sec"<<endl;
		oracle_TimeFile.close(); 
	}		

	/*#ifdef PRUNE_STATS
		cout<<"numBound: "<<stat.numBound<<endl;
		cout<<"numExact: "<<stat.numExact<<endl;
	#endif

	#ifdef COMPONENT_CLOCK
		cout<<"BC Time (on average): "<<stat.BC_Time/(double)stat.numBound<<endl;
		cout<<"EC Time (on average): "<<stat.EC_Time/(double)stat.numExact<<endl;
	#endif*/

	/*if((method>=13 && method<=19) || method==101)
	{
		cout<<"Statisitcs of each level:"<<endl;

		cout<<(double)stat.incr_Level_Counter[(int)stat.incr_Level_Counter.size()-2]/qNum<<" ";
		cout<<(double)stat.incr_Level_Counter[(int)stat.incr_Level_Counter.size()-1]/qNum<<" ";

		for(int l=0;l<(int)stat.incr_Level_Counter.size()-2;l++)
			cout<<(double)stat.incr_Level_Counter[l]/qNum<<" ";
		cout<<endl;
	}*/
}
