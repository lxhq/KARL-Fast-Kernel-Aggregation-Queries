#include "Incremental_Classification.h"

struct group_Entry
{
	int groupNum;
	double LB;
	double UB;
};

//This is the maximum heap
struct compare_Group_Priority
{
	bool operator()(group_Entry& e1, group_Entry& e2)
	{
		return (e1.UB-e1.LB)<(e2.UB-e2.LB);
	}
};

typedef priority_queue<group_Entry,vector<group_Entry>,compare_Group_Priority> group_Heap;

double LB_T(double*q,group_Info& g_Info,double& gamma_sum,int dim,SVM_stat& stat)
{
	double ip=0;
	for(int d=0;d<dim;d++)
		ip+=q[d]*g_Info.a_G[d];

	gamma_sum=stat.gammaValue*(g_Info.w_G*stat.qSquareNorm-2*ip+g_Info.S_G);

	return g_Info.w_G*exp(-(gamma_sum/g_Info.w_G));
}

double UB_C(double*q,group_Info& g_Info,double& gamma_sum,int dim,SVM_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;
	double temp_obt_dist;

	//l_tri and u_tri
	ub=u_tri(q,g_Info.center,dim,g_Info.radius,temp_obt_dist);

	if(temp_obt_dist<g_Info.radius)
		lb=0;
	else
		lb=temp_obt_dist-g_Info.radius;

	l2=lb*lb;
	u2=ub*ub;

	//l2=u2, we cannot use mx+c to act as upper bound computationally
	if(u2-l2<epsilon)
		return g_Info.w_G*exp(-stat.gammaValue*l2);

	exp_l2=exp(-stat.gammaValue*l2);
	exp_u2=exp(-stat.gammaValue*u2);

	//compute m and c
	m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*g_Info.w_G);
}

void incremental_Classification(double*q,double**P,group_Info*groupArray,double*alphaArray,int chunkSize,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q)
{
	double curLB;
	double curUB;
	double tempExactValue;
	int num_of_Incr;
	double gamma_sum;

	//num_of_Incr=(int)ceil((double)stat.total_sv/chunkSize);
	num_of_Incr=(int)floor((double)stat.total_sv/chunkSize);

	double*tempLB_Array=new double[num_of_Incr];
	double*tempUB_Array=new double[num_of_Incr];

	curLB=LB_T(q,groupArray[num_of_Incr],gamma_sum,dim,stat);

	if(curLB>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr]++;
		return;
	}

	curUB=UB_C(q,groupArray[num_of_Incr],gamma_sum,dim,stat);

	if(curUB<=stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr]++;
		return;
	}

	//Obtain (LB/UB) bound values for each partition
	//***************************************************************//
	curLB=0;
	curUB=0;
	for(int i=0;i<num_of_Incr;i++)
	{
		tempLB_Array[i]=LB_T(q,groupArray[i],gamma_sum,dim,stat);
		tempUB_Array[i]=UB_C(q,groupArray[i],gamma_sum,dim,stat);
		curLB=curLB+tempLB_Array[i];
		curUB=curUB+tempUB_Array[i];
	}

	if(curLB>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr+1]++;
		return;
	}

	if(curUB<=stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr+1]++;
		return;
	}
	//***************************************************************//

	//Incremental Part
	for(int i=0;i<num_of_Incr-1;i++)
	{
		curLB-=tempLB_Array[i];
		curUB-=tempUB_Array[i];

		tempExactValue=callSVM_Refine(svm_q,model,i*chunkSize,(i+1)*chunkSize);

		curLB+=tempExactValue;

		if(curLB>=stat.rho)
		{
			stat.class_resultVector.push_back(1);
			stat.pruneCount++;
			stat.incr_Level_Counter[i]++;
			return;
		}

		curUB+=tempExactValue;

		if(curUB<=stat.rho)
		{
			stat.class_resultVector.push_back(-1);
			stat.pruneCount++;
			stat.incr_Level_Counter[i]++;
			return;
		}
	}

	curLB-=tempLB_Array[num_of_Incr-1];
	curUB-=tempUB_Array[num_of_Incr-1];

	int lastChunkInit=(num_of_Incr-1)*chunkSize;
	tempExactValue=callSVM_Refine(svm_q,model,lastChunkInit,stat.total_sv);

	curLB+=tempExactValue;
	curUB+=tempExactValue;

	if(curLB>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.incr_Level_Counter[num_of_Incr-1]++;
		return;
	}
	if(curUB<=stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.incr_Level_Counter[num_of_Incr-1]++;
		return;
	}
}

void incremental_Classification_multiStep(double*q,double**P,group_Info*groupArray,double*alphaArray,int chunkSize,int dim,SVM_stat& stat,svm_model*model,svm_node*svm_q)
{
	double curLB;
	double curUB;
	double tempExactValue;
	int num_of_Incr;
	double gamma_sum;
	group_Heap g_Heap;
	group_Entry g_Entry;

	//statistical information
	int level;

	//num_of_Incr=(int)ceil((double)stat.total_sv/chunkSize);
	num_of_Incr=(int)floor((double)stat.total_sv/chunkSize);

	curLB=LB_T(q,groupArray[num_of_Incr],gamma_sum,dim,stat);

	if(curLB>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr]++;
		return;
	}

	curUB=UB_C(q,groupArray[num_of_Incr],gamma_sum,dim,stat);

	if(curUB<=stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr]++;
		return;
	}

	//Obtain (LB/UB) bound values for each partition
	//***************************************************************//
	curLB=0;
	curUB=0;
	for(int i=0;i<num_of_Incr;i++)
	{
		g_Entry.LB=LB_T(q,groupArray[i],gamma_sum,dim,stat);
		g_Entry.UB=UB_C(q,groupArray[i],gamma_sum,dim,stat);
		g_Entry.groupNum=i;

		g_Heap.push(g_Entry);

		curLB=curLB+g_Entry.LB;
		curUB=curUB+g_Entry.UB;
	}

	if(curLB>=stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr+1]++;
		return;
	}

	if(curUB<=stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		stat.incr_Level_Counter[num_of_Incr+1]++;
		return;
	}
	//***************************************************************//

	level=0;
	while(g_Heap.size()>0)
	{
		g_Entry=g_Heap.top();
		g_Heap.pop();

		if(g_Entry.groupNum!=num_of_Incr-1)
			tempExactValue=callSVM_Refine(svm_q,model,g_Entry.groupNum*chunkSize,(g_Entry.groupNum+1)*chunkSize);
		else
			tempExactValue=callSVM_Refine(svm_q,model,(num_of_Incr-1)*chunkSize,stat.total_sv);

		curLB=curLB-g_Entry.LB+tempExactValue;

		if(curLB>=stat.rho)
		{
			stat.class_resultVector.push_back(1);
			stat.pruneCount++;
			stat.incr_Level_Counter[level]++;
			return;
		}

		curUB=curUB-g_Entry.UB+tempExactValue;

		if(curUB<=stat.rho)
		{
			stat.class_resultVector.push_back(-1);
			stat.pruneCount++;
			stat.incr_Level_Counter[level]++;
			return;
		}

		level++;
	}
}

void update_group(double**matrix,double*sub_alphaArray,group_Info*groupArray,int g,int max_size,int dim)
{
	groupArray[g].a_G=new double[dim];
	groupArray[g].center=new double[dim];
	groupArray[g].S_G=0;
	groupArray[g].w_G=0;
	groupArray[g].radius=0;

	double temp_Radius;
	double sq_norm;

	for(int d=0;d<dim;d++)
		groupArray[g].a_G[d]=0;
	
	for(int i=0;i<max_size;i++)
	{
		sq_norm=0;
		groupArray[g].w_G+=sub_alphaArray[i];
		for(int d=0;d<dim;d++)
		{
			sq_norm+=matrix[i][d]*matrix[i][d];
			groupArray[g].a_G[d]+=sub_alphaArray[i]*matrix[i][d];
		}
		groupArray[g].S_G+=sub_alphaArray[i]*sq_norm;
	}

	for(int d=0;d<dim;d++)
		groupArray[g].center[d]=groupArray[g].a_G[d]/groupArray[g].w_G;

	for(int i=0;i<max_size;i++)
	{
		temp_Radius=0;
		for(int d=0;d<dim;d++)
			temp_Radius+=(matrix[i][d]-groupArray[g].center[d])*(matrix[i][d]-groupArray[g].center[d]);

		if(temp_Radius>groupArray[g].radius)
			groupArray[g].radius=temp_Radius;
	}
	groupArray[g].radius=sqrt(groupArray[g].radius);
}

void create_groupArray(double**P,double*alphaArray,group_Info*& groupArray,int chunkSize,int dim,SVM_stat& stat,bool initPruneCounter)
{
	//int num_of_incr=(int)ceil((double)stat.total_sv/chunkSize);
	int num_of_incr=(int)floor((double)stat.total_sv/chunkSize);
	int max_size;

	groupArray=new group_Info[num_of_incr+1];

	update_group(&P[0],&alphaArray[0],groupArray,num_of_incr,stat.total_sv,dim);

	for(int g=num_of_incr-1;g>=0;g--)
	{
		if(g==num_of_incr-1)
			max_size=stat.total_sv-g*chunkSize;
		else
			max_size=chunkSize;

		update_group(&P[g*chunkSize],&alphaArray[g*chunkSize],groupArray,g,max_size,dim);
	}

	//add for detecting the statistics of each level
	//***********************************************//
	if(initPruneCounter==true)
	{
		stat.incr_Level_Counter.push_back(0);
		stat.incr_Level_Counter.push_back(0);
		for(int g=num_of_incr-1;g>=0;g--)
			stat.incr_Level_Counter.push_back(0);
	}
	//***********************************************//
}