#include "Sequential_Bound.h"

void pre_Compute_sequential(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double*& a_G,double& S_G,double*& center,double& radius,SVM_stat& stat)
{
	double squareNorm;
	double tempRadius;
	//initialization
	sum_alpha=0;
	a_G=new double[dim];
	center=new double[dim];
	for(int d=0;d<dim;d++)
	{
		a_G[d]=0;
		center[d]=0;
	}
	S_G=0;
	radius=-inf;

	for(int i=0;i<stat.total_sv;i++)
	{
		squareNorm=0;
		for(int d=0;d<dim;d++)
		{
			a_G[d]+=alphaArray[i]*dataMatrix[i][d];
			squareNorm+=dataMatrix[i][d]*dataMatrix[i][d];
		}
		S_G+=alphaArray[i]*squareNorm;
		sum_alpha+=alphaArray[i];
	}

	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sum_alpha;

	for(int i=0;i<stat.total_sv;i++)
	{
		tempRadius=0;
		for(int d=0;d<dim;d++)
			tempRadius+=(dataMatrix[i][d]-center[d])*(dataMatrix[i][d]-center[d]);

		if(tempRadius>radius)
			radius=tempRadius;
	}

	radius=sqrt(radius);
}

void sequential_bound(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,int method)
{
	double LB,UB;
	double ell,u;
	double m,c;
	double obt_dist;
	double sum_of_sqEuclid;
	double q_sq_norm=0;
	double ip=0;
	double exp_ell2,exp_u2;

	for(int d=0;d<dim;d++)
	{
		q_sq_norm+=q[d]*q[d];
		ip+=q[d]*a_G[d];
	}

	sum_of_sqEuclid=sum_alpha*q_sq_norm-2*ip+S_G;

	//Compute lower bound
	//LB=stat.gammaValue/sum_alpha*sum_of_sqEuclid;
	LB=sum_alpha*exp(-stat.gammaValue/sum_alpha*sum_of_sqEuclid);

	if(LB>stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		return;
	}

	//Compute upper bound
	u=u_tri(q,center,dim,radius,obt_dist);
	ell=obt_dist-radius;
	if(ell<0)
		ell=0;

	exp_ell2=exp(-stat.gammaValue*ell*ell);
	exp_u2=exp(-stat.gammaValue*u*u);
	m=(exp_u2-exp_ell2)/(stat.gammaValue*(u*u-ell*ell));
	c=(u*u*exp_ell2-ell*ell*exp_u2)/(u*u-ell*ell);

	UB=m*stat.gammaValue*sum_of_sqEuclid+c*sum_alpha;
	if(UB<stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		return;
	}

	if(method==9)
		prediction(model,svm_q,stat);
	if(method==10)
		SS_iter(q,dataMatrix,alphaArray,dim,stat);
}

void pre_Compute_sequential_MBR(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double**& boundary,SVM_stat& stat)
{
	sum_alpha=0;
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
		for(int i=0;i<(int)stat.total_sv;i++)
		{
			if(dataMatrix[i][d]<boundary[d][0])
				boundary[d][0]=dataMatrix[i][d];

			if(dataMatrix[i][d]>boundary[d][1])
				boundary[d][1]=dataMatrix[i][d];
		}
	}

	for(int i=0;i<(int)stat.total_sv;i++)
		sum_alpha+=alphaArray[i];
}

void sequential_bound_MBR(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double**boundary,SVM_stat& stat,int method)
{
	double LB,UB;
	double ell,u;

	ell=ell_MBR(q,boundary,dim);
	u=u_MBR(q,boundary,dim);

	LB=sum_alpha*exp(-stat.gammaValue*u*u);
	UB=sum_alpha*exp(-stat.gammaValue*ell*ell);

	if(LB>stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		return;
	}
	if(UB<stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		return;
	}

	if(method==30)
		prediction(model,svm_q,stat);
	if(method==31)
		SS_iter(q,dataMatrix,alphaArray,dim,stat);
}

void pre_Compute_sequential_Delta(double**dataMatrix,double*alphaArray,int dim,double& sum_alpha,double*& center,double& radius,SVM_stat& stat)
{
	double tempRadius;
	//initialization
	sum_alpha=0;
	center=new double[dim];
	for(int d=0;d<dim;d++)
		center[d]=0;

	for(int i=0;i<stat.total_sv;i++)
	{
		for(int d=0;d<dim;d++)
			center[d]+=alphaArray[i]*dataMatrix[i][d];
		sum_alpha+=alphaArray[i];
	}

	for(int d=0;d<dim;d++)
		center[d]=center[d]/sum_alpha;

	radius=-inf;
	for(int i=0;i<stat.total_sv;i++)
	{
		tempRadius=0;
		for(int d=0;d<dim;d++)
			tempRadius+=(dataMatrix[i][d]-center[d])*(dataMatrix[i][d]-center[d]);

		if(tempRadius>radius)
			radius=tempRadius;
	}

	radius=sqrt(radius);
}

void sequential_bound_Delta(double*q,double**dataMatrix,double*alphaArray,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha,double*center,double radius,SVM_stat& stat,int method)
{
	double LB,UB;
	double ell,u;
	double temp_obt_dist;

	//l_tri and u_tri
	u=u_tri(q,center,dim,radius,temp_obt_dist);
	if(temp_obt_dist<radius)
		ell=0;
	else
		ell=temp_obt_dist-radius;

	LB=sum_alpha*exp(-stat.gammaValue*u*u);
	UB=sum_alpha*exp(-stat.gammaValue*ell*ell);

	if(LB>stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		return;
	}
	if(UB<stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		return;
	}

	if(method==32)
		prediction(model,svm_q,stat);
	if(method==33)
		SS_iter(q,dataMatrix,alphaArray,dim,stat);
}

void sequential_bound_sparse(double*q,double*alphaArray,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,double*p_sqNormArray,double**lu_Table)
{
	double LB,UB;
	double ell,u;
	double m,c;
	double obt_dist;
	double sum_of_sqEuclid;
	double q_sq_norm=0;
	double ip=0;
	double exp_ell2,exp_u2;

	for(int d=0;d<dim;d++)
	{
		q_sq_norm+=q[d]*q[d];
		ip+=q[d]*a_G[d];
	}

	sum_of_sqEuclid=sum_alpha*q_sq_norm-2*ip+S_G;

	//Compute lower bound
	//LB=stat.gammaValue/sum_alpha*sum_of_sqEuclid;
	LB=sum_alpha*exp(-stat.gammaValue/sum_alpha*sum_of_sqEuclid);

	if(LB>stat.rho)
	{
		stat.class_resultVector.push_back(1);
		stat.pruneCount++;
		return;
	}

	//Compute upper bound
	u=u_tri(q,center,dim,radius,obt_dist);
	ell=obt_dist-radius;
	if(ell<0)
		ell=0;

	exp_ell2=exp(-stat.gammaValue*ell*ell);
	exp_u2=exp(-stat.gammaValue*u*u);
	m=(exp_u2-exp_ell2)/(stat.gammaValue*(u*u-ell*ell));
	c=(u*u*exp_ell2-ell*ell*exp_u2)/(u*u-ell*ell);

	UB=m*stat.gammaValue*sum_of_sqEuclid+c*sum_alpha;
	if(UB<stat.rho)
	{
		stat.class_resultVector.push_back(-1);
		stat.pruneCount++;
		return;
	}

	SS_sparse_iter(q,alphaArray,dim,stat,inv_Index,sq_EuclidVector,q_sq_norm,p_sqNormArray,lu_Table);
}