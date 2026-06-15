#include "SS.h"

void SS_iter(double*q,double**dataMatrix,int dim,KDE_stat& stat)
{
	double incr_Value=0;
	double sq_Value;

	for(int i=0;i<stat.n;i++)
	{
		sq_Value=0;
		for(int d=0;d<dim;d++)
			sq_Value+=(q[d]-dataMatrix[i][d])*(q[d]-dataMatrix[i][d]);

		incr_Value=incr_Value+exp(-sq_Value);
	}

	stat.resultValueVector.push_back(incr_Value);
}

void pre_Compute_sequential_MBR(double**dataMatrix,int dim,double**& boundary,KDE_stat& stat)
{
	boundary=new double*[dim];
	for(int d=0;d<dim;d++)
		boundary[d]=new double[2];

	for(int d=0;d<dim;d++)
	{
		boundary[d][0]=inf;
		boundary[d][1]=-inf;
	}

	for(int d=0;d<dim;d++)
	{
		for(int i=0;i<(int)stat.n;i++)
		{
			if(dataMatrix[i][d]<boundary[d][0])
				boundary[d][0]=dataMatrix[i][d];

			if(dataMatrix[i][d]>boundary[d][1])
				boundary[d][1]=dataMatrix[i][d];
		}
	}
}

void SS_MBR(double*q,double**dataMatrix,int dim,double**boundary,KDE_stat& stat)
{
	double LB,UB;
	double ell,u;
	double validate_Value;

	ell=ell_MBR(q,boundary,dim);
	u=u_MBR(q,boundary,dim);

	LB=stat.n*exp(-u*u);
	UB=stat.n*exp(-ell*ell);

	if(validate_best(LB,UB,stat.rel_error,validate_Value)==true)
		stat.resultValueVector.push_back(validate_Value);
	else
		SS_iter(q,dataMatrix,dim,stat);
}

void pre_Compute_sequential_Delta(double**dataMatrix,int dim,double*& center,double& radius,KDE_stat& stat)
{
	double tempRadius;
	//initialization
	center=new double[dim];
	for(int d=0;d<dim;d++)
		center[d]=0;

	for(int i=0;i<stat.n;i++)
	{
		for(int d=0;d<dim;d++)
			center[d]+=dataMatrix[i][d];
	}

	for(int d=0;d<dim;d++)
		center[d]=center[d]/((double)stat.n);

	radius=-inf;
	for(int i=0;i<stat.n;i++)
	{
		tempRadius=0;
		for(int d=0;d<dim;d++)
			tempRadius+=(dataMatrix[i][d]-center[d])*(dataMatrix[i][d]-center[d]);

		if(tempRadius>radius)
			radius=tempRadius;
	}

	radius=sqrt(radius);
}

void SS_Delta(double*q,double**dataMatrix,int dim,double*center,double radius,KDE_stat& stat)
{
	double LB,UB;
	double ell,u;
	double temp_obt_dist;
	double validate_Value;

	//l_tri and u_tri
	u=u_tri(q,center,dim,radius,temp_obt_dist);
	if(temp_obt_dist<radius)
		ell=0;
	else
		ell=temp_obt_dist-radius;

	LB=stat.n*exp(-u*u);
	UB=stat.n*exp(-ell*ell);

	if(validate_best(LB,UB,stat.rel_error,validate_Value)==true)
		stat.resultValueVector.push_back(validate_Value);
	else
		SS_iter(q,dataMatrix,dim,stat);
}

void pre_Compute_sequential(double**dataMatrix,int dim,double**& boundary,double*& a_G,double& S_G,KDE_stat& stat)
{
	double squareNorm;

	//initialization
	a_G=new double[dim];
	for(int d=0;d<dim;d++)
		a_G[d]=0;
	S_G=0;

	for(int i=0;i<stat.n;i++)
	{
		squareNorm=0;
		for(int d=0;d<dim;d++)
		{
			a_G[d]+=dataMatrix[i][d];
			squareNorm+=dataMatrix[i][d]*dataMatrix[i][d];
		}
		S_G+=squareNorm;
	}

	pre_Compute_sequential_MBR(dataMatrix,dim,boundary,stat);
}

void SS_linear(double*q,double**dataMatrix,int dim,double**& boundary,double*a_G,double S_G,KDE_stat& stat)
{
	double LB,UB;
	double ell,u;
	double m,c;
	double sum_of_sqEuclid;
	double q_sq_norm=0;
	double ip=0;
	double exp_ell2,exp_u2;
	double validate_Value;

	for(int d=0;d<dim;d++)
	{
		q_sq_norm+=q[d]*q[d];
		ip+=q[d]*a_G[d];
	}

	sum_of_sqEuclid=stat.n*q_sq_norm-2*ip+S_G;

	//Compute lower bound
	LB=stat.n*exp(-sum_of_sqEuclid/((double)stat.n));

	//Compute upper bound
	ell=ell_MBR(q,boundary,dim);
	u=u_MBR(q,boundary,dim);

	exp_ell2=exp(-ell*ell);
	exp_u2=exp(-u*u);
	m=(exp_u2-exp_ell2)/(u*u-ell*ell);
	c=(u*u*exp_ell2-ell*ell*exp_u2)/(u*u-ell*ell);

	UB=m*sum_of_sqEuclid+c*stat.n;

	if(validate_best(LB,UB,stat.rel_error,validate_Value)==true)
		stat.resultValueVector.push_back(validate_Value);
	else
		SS_iter(q,dataMatrix,dim,stat);
}