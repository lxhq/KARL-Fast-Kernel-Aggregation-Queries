#include "exp_LookUp.h"

void lookUpTable_Creation(double**& lu_Table)
{
	//Initialization
	lu_Table=new double*[TableSize];
	for(int t=0;t<TableSize;t++)
		lu_Table[t]=new double[2];

	double t_double;
	for(int t=0;t<TableSize;t++)
	{
		t_double=((double)t)*0.000001;
		lu_Table[t][0]=exp(-(t_double+0.000001));//LB
		lu_Table[t][1]=exp(-t_double);//UB
	}
}

double lookUpLB(double x,double**lu_Table)
{
	int index=(int)(x*1000000);
	//added to prevent the negative index
	if(index<0)
		return -1;

	if(index>=TableSize)
		return 0;
	else
		return lu_Table[index][0];
}

double lookUpUB(double x,double**lu_Table)
{
	int index=(int)(x*1000000);
	if(index>=TableSize)
		return lu_Table[TableSize-1][1];
	else
		return lu_Table[index][1];
}

void lookUp_sim_Bounds(double x,double& l,double& u,double**lu_Table)
{
	int index=(int)(x*1000000);
	if(index>TableSize)
	{
		l=0;
		u=lu_Table[TableSize-1][1];
	}
	else
	{
		l=lu_Table[index][0];
		u=lu_Table[index][1];
	}
}