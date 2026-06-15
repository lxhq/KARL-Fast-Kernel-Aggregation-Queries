#include "Oracle.h"

//Online Oracle
void clear_groupArray(group_Info*groupArray,SVM_stat& stat,int chunkSize)
{
	int num_of_incr=(int)floor((double)stat.total_sv/chunkSize);
	for(int i=0;i<num_of_incr+1;i++)
	{
		delete[] groupArray[i].a_G;
		delete[] groupArray[i].center;
	}
	delete[] groupArray;
}

double Oracle_Algorithm(double*q,double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q)
{
	group_Info*groupArray;

	//clock_t start_s;
	//clock_t stop_s;

	int max_ChunkSize=(int)((double)stat.total_sv/2.0);
	double*time_counter=new double[10];
	int chunkSize=0;

	for(int chunkNum=0;chunkNum<10;chunkNum++)
	{
		chunkSize=(int)((((double)chunkNum+1.0)/10.0)*max_ChunkSize);

		create_groupArray(P,alphaArray,groupArray,chunkSize,dim,stat,true);

		#ifdef C_PLUSPLUS11_CLOCK
			auto start_s=chrono::high_resolution_clock::now();
		#else
			clock_t start_s=clock();
		#endif
		
		incremental_Classification_multiStep(q,P,groupArray,alphaArray,chunkSize,dim,stat,model,svm_q);
		
		#ifdef C_PLUSPLUS11_CLOCK
			auto stop_s=chrono::high_resolution_clock::now();
			time_counter[chunkNum]=chrono::duration_cast<chrono::nanoseconds>(stop_s-start_s).count();
		#else
			clock_t stop_s=clock();
			time_counter[chunkNum]=(double)(stop_s-start_s);
		#endif
		
		clear_groupArray(groupArray,stat,chunkSize);
		stat.incr_Level_Counter.clear();
	}

	double bestTime=inf;
	//for(int c=2;c<=max_ChunkSize;c++)
	for(int chunkNum=0;chunkNum<10;chunkNum++)
	{
		if(time_counter[chunkNum]<bestTime)
			#ifdef C_PLUSPLUS11_CLOCK
				bestTime=time_counter[chunkNum];
			#else
				bestTime=time_counter[chunkNum]/CLOCKS_PER_SEC;
			#endif
	}

	delete[] time_counter;

	return bestTime;
}

double Oracle_Weak_Algorithm(double*q,double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q)
{
	static int Counter=0;
	group_Info*groupArray;

	//clock_t start_s;
	//clock_t stop_s;
	double bestTime;

	//int max_ChunkSize=(int)((double)stat.total_sv/2.0);
	int chunkSize=(int)(1.0/5.0*stat.total_sv);
	if(Counter>0)
		create_groupArray(P,alphaArray,groupArray,chunkSize,dim,stat,false);
	else
	{
		create_groupArray(P,alphaArray,groupArray,chunkSize,dim,stat,true);
		Counter++;
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_s=clock();
	#endif
	
	incremental_Classification_multiStep(q,P,groupArray,alphaArray,chunkSize,dim,stat,model,svm_q);
	
	#ifdef C_PLUSPLUS11_CLOCK
		auto stop_s=chrono::high_resolution_clock::now();
		bestTime=chrono::duration_cast<chrono::nanoseconds>(stop_s-start_s).count();
	#else
		clock_t stop_s=clock();
		bestTime=(double)(stop_s-start_s)/CLOCKS_PER_SEC;
	#endif

	return bestTime;
}

void convert_q_to_SVMFormat(double**queryMatrix,svm_node*& svm_q,int select_q,int dim)
{
	int position=0;
	svm_q=new svm_node[dim];

	for(int d=0;d<dim;d++)
	{
		if(fabs(queryMatrix[select_q][d])>epsilon)
		{
			svm_q[position].index=d+1;
			svm_q[position].value=queryMatrix[select_q][d];
			position++;
		}
	}
	svm_q[position].index=-1;
}

double runOracle(double*q,double**P,svm_node*svm_q,svm_model*model,double*alphaArray,int dim,SVM_stat& stat,svm_node*x_space,int method)
{
	double best_time;
	vector<orderEntry> orderList;

	//Reordering
	oracle_Reorder(q,P,orderList,dim,stat);
	reordering(P,alphaArray,orderList,stat,dim);
	createModel_inMemory(P,alphaArray,dim,stat,model,x_space);

	if(method==100)
		best_time=Oracle_Algorithm(q,P,alphaArray,dim,stat,model,svm_q);
	else //method=101
		best_time=Oracle_Weak_Algorithm(q,P,alphaArray,dim,stat,model,svm_q);

	orderList.clear();
	return best_time;
}

void outputOracle_QueryTime(char*Oracle_TimeFileName,double bestTime)
{
	fstream oracle_TimeFile;

	oracle_TimeFile.open(Oracle_TimeFileName,ios::out|ios::app);

	if(oracle_TimeFile.is_open()==false)
	{
		cout<<"Cannot Open oracle_TimeFile!"<<endl;
		exit(1);
	}

	oracle_TimeFile<<bestTime<<endl;
	oracle_TimeFile.close();
}