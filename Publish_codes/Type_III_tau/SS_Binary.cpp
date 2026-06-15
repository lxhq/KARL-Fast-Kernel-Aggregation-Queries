#include "SS_Binary.h"

void SS_iter(double*q,double**dataMatrix,double*outputArray,int dim,SVM_stat& stat)
{
	double incr_Value=0;
	double sq_Value;

	for(int i=0;i<stat.total_sv;i++)
	{
		sq_Value=0;
		for(int d=0;d<dim;d++)
			sq_Value+=(q[d]-dataMatrix[i][d])*(q[d]-dataMatrix[i][d]);

		incr_Value=incr_Value+outputArray[i]*exp(-stat.gammaValue*sq_Value);
	}

	#ifdef EXACT_VALUE_STATS
		stat.exactValueVector.push_back(incr_Value);
	#endif

	if(incr_Value>stat.rho)
		stat.class_resultVector.push_back(1);
	else
		stat.class_resultVector.push_back(-1);
}


void pre_Compute_Bplus(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double*& a_G_pos,double*& a_G_neg,double& S_G_pos,double& S_G_neg,double*& center_pos,double*& center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat)
{
	double squareNorm;
	double tempRadius;

	//Initialization
	sum_alpha_pos=0;
	sum_alpha_neg=0;
	S_G_pos=0;
	S_G_neg=0;
	a_G_pos=new double[dim];
	a_G_neg=new double[dim];
	center_pos=new double[dim];
	center_neg=new double[dim];
	radius_pos=-inf;
	radius_neg=-inf;

	for(int d=0;d<dim;d++)
	{
		a_G_pos[d]=0;
		a_G_neg[d]=0;
	}

	//update a_G,S_G and sum_alpha (w_G)
	for(int i=0;i<stat.total_sv;i++)
	{
		squareNorm=0;
		for(int d=0;d<dim;d++)
			squareNorm+=dataMatrix[i][d]*dataMatrix[i][d];

		if(outputArray[i]>0)
		{
			for(int d=0;d<dim;d++)
				a_G_pos[d]+=outputArray[i]*dataMatrix[i][d];

			S_G_pos+=outputArray[i]*squareNorm;
			sum_alpha_pos+=outputArray[i];
		}
		else
		{
			for(int d=0;d<dim;d++)
				a_G_neg[d]+=fabs(outputArray[i])*dataMatrix[i][d];

			S_G_neg+=fabs(outputArray[i])*squareNorm;
			sum_alpha_neg+=fabs(outputArray[i]);
		}
	}

	for(int d=0;d<dim;d++)
	{
		center_pos[d]=a_G_pos[d]/sum_alpha_pos;
		center_neg[d]=a_G_neg[d]/sum_alpha_neg;
	}

	//update radius
	for(int i=0;i<stat.total_sv;i++)
	{
		if(outputArray[i]>0)
		{
			tempRadius=euclid_dist(dataMatrix[i],center_pos,dim);
			if(tempRadius>radius_pos)
				radius_pos=tempRadius;
		}
		else
		{
			tempRadius=euclid_dist(dataMatrix[i],center_neg,dim);
			if(tempRadius>radius_neg)
				radius_neg=tempRadius;
		}
	}

}

void LibSVM_Bplus(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double*& a_G_pos,double*& a_G_neg,double& S_G_pos,double& S_G_neg,double*center_pos,double*center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat)
{
	double LB,UB;
	double ell,u;
	double m,c;
	double obt_dist;
	double sum_of_sqEuclid;
	double q_sq_norm=0;
	double ip_neg=0;
	double ip_pos=0;

	double exp_ell2,exp_u2;

	double LB_pos, UB_pos;
	double LB_neg, UB_neg;

	for(int d=0;d<dim;d++)
	{
		q_sq_norm+=q[d]*q[d];
		ip_pos+=q[d]*a_G_pos[d];
		ip_neg+=q[d]*a_G_neg[d];
	}

	//obtain LB_pos and UB_pos 
	//**************************************************************************//
	sum_of_sqEuclid=sum_alpha_pos*q_sq_norm-2*ip_pos+S_G_pos;
	LB_pos=sum_alpha_pos*exp(-stat.gammaValue/sum_alpha_pos*sum_of_sqEuclid);

	u=u_tri(q,center_pos,dim,radius_pos,obt_dist);
	ell=obt_dist-radius_pos;
	if(ell<0)
		ell=0;

	exp_ell2=exp(-stat.gammaValue*ell*ell);
	exp_u2=exp(-stat.gammaValue*u*u);
	m=(exp_u2-exp_ell2)/(stat.gammaValue*(u*u-ell*ell));
	c=(u*u*exp_ell2-ell*ell*exp_u2)/(u*u-ell*ell);

	UB_pos=m*stat.gammaValue*sum_of_sqEuclid+c*sum_alpha_pos;
	//**************************************************************************//


	//obtain LB_neg and UB_neg
	//**************************************************************************//
	sum_of_sqEuclid=sum_alpha_neg*q_sq_norm-2*ip_neg+S_G_neg;
	LB_neg=sum_alpha_neg*exp(-stat.gammaValue/sum_alpha_neg*sum_of_sqEuclid);

	u=u_tri(q,center_neg,dim,radius_neg,obt_dist);
	ell=obt_dist-radius_neg;
	if(ell<0)
		ell=0;

	exp_ell2=exp(-stat.gammaValue*ell*ell);
	exp_u2=exp(-stat.gammaValue*u*u);
	m=(exp_u2-exp_ell2)/(stat.gammaValue*(u*u-ell*ell));
	c=(u*u*exp_ell2-ell*ell*exp_u2)/(u*u-ell*ell);

	UB_neg=m*stat.gammaValue*sum_of_sqEuclid+c*sum_alpha_neg;
	//**************************************************************************//

	LB=LB_pos-UB_neg;
	UB=UB_pos-LB_neg;

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

	prediction(model,svm_q,stat);
}

//LibSVM_{MBR}
void pre_Compute_MBR(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double**& boundary_pos,double**& boundary_neg,SVM_stat& stat)
{
	sum_alpha_pos=0;
	sum_alpha_neg=0;
	boundary_pos=new double*[dim];
	boundary_neg=new double*[dim];
	for(int d=0;d<dim;d++)
	{
		boundary_pos[d]=new double[2];
		boundary_neg[d]=new double[2];
	}

	for(int d=0;d<dim;d++)
	{
		boundary_pos[d][0]=inf;
		boundary_pos[d][1]=-inf;

		boundary_neg[d][0]=inf;
		boundary_neg[d][1]=-inf;
	}

	
	for(int i=0;i<(int)stat.total_sv;i++)
	{
		if(outputArray[i]>0)
		{
			sum_alpha_pos+=outputArray[i];
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[i][d]<boundary_pos[d][0])
					boundary_pos[d][0]=dataMatrix[i][d];

				if(dataMatrix[i][d]>boundary_pos[d][1])
					boundary_pos[d][1]=dataMatrix[i][d];
			}
		}
		else
		{
			sum_alpha_neg+=fabs(outputArray[i]);
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[i][d]<boundary_neg[d][0])
					boundary_neg[d][0]=dataMatrix[i][d];

				if(dataMatrix[i][d]>boundary_neg[d][1])
					boundary_neg[d][1]=dataMatrix[i][d];
			}
		}
	}

}

void LibSVM_MBR(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double**& boundary_pos,double**& boundary_neg,SVM_stat& stat)
{
	double lb_neg;
	double ub_pos;
	double L_pos;
	double U_neg;

	double lb_pos;
	double ub_neg;
	double L_neg;
	double U_pos;
	double LB;
	double UB;

	if(sum_alpha_pos<epsilon)
		ub_pos=0;
	else
		ub_pos=u_MBR(q,boundary_pos,dim); //G_{+}	

	if(sum_alpha_neg<epsilon)
		lb_neg=0;
	else
		lb_neg=ell_MBR(q,boundary_neg,dim); //G_{-}


	if(sum_alpha_pos<epsilon)
		lb_pos=0;
	else
		lb_pos=ell_MBR(q,boundary_pos,dim); //G_{+}

	if(sum_alpha_neg<epsilon)
		ub_neg=0;
	else
		ub_neg=u_MBR(q,boundary_neg,dim); //G_{-}

	L_pos=sum_alpha_pos*exp(-stat.gammaValue*ub_pos*ub_pos); //LB_{G_{+}}
	U_neg=sum_alpha_neg*exp(-stat.gammaValue*lb_neg*lb_neg); //UB_{G_{-}}

	L_neg=sum_alpha_neg*exp(-stat.gammaValue*ub_neg*ub_neg); //LB_{G_{-}}
	U_pos=sum_alpha_pos*exp(-stat.gammaValue*lb_pos*lb_pos); //UB_{G_{+}}

	LB=L_pos-U_neg;
	UB=U_pos-L_neg;

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

	prediction(model,svm_q,stat);
}

void pre_Compute_Delta(double**dataMatrix,double*outputArray,int dim,double& sum_alpha_pos,double& sum_alpha_neg,double*& center_pos,double*& center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat)
{
	double dist;
	center_pos=new double[dim];
	center_neg=new double[dim];
	for(int d=0;d<dim;d++)
	{
		center_pos[d]=0;
		center_neg[d]=0;
	}
	sum_alpha_pos=0;
	sum_alpha_neg=0;

	for(int i=0;i<(int)stat.total_sv;i++)
	{
		for(int d=0;d<dim;d++)
		{
			if(outputArray[i]>=0)
				center_pos[d]=center_pos[d]+outputArray[i]*dataMatrix[i][d];
			else
				center_neg[d]=center_neg[d]+fabs(outputArray[i])*dataMatrix[i][d];
		}

		if(outputArray[i]>=0)
			sum_alpha_pos+=outputArray[i];
		else
			sum_alpha_neg+=fabs(outputArray[i]);
	}

	for(int d=0;d<dim;d++)
	{
		center_pos[d]=center_pos[d]/sum_alpha_pos;
		center_neg[d]=center_neg[d]/sum_alpha_neg;
	}

	//update radius
	radius_pos=0;
	radius_neg=0;
	for(int i=0;i<(int)stat.total_sv;i++)
	{
		if(outputArray[i]>=0)
		{
			dist=euclid_dist(center_pos,dataMatrix[i],dim);
			if(dist>radius_pos)
				radius_pos=dist;
		}
		else
		{
			dist=euclid_dist(center_neg,dataMatrix[i],dim);
			if(dist>radius_neg)
				radius_neg=dist;
		}
	}
}

void LibSVM_Delta(double*q,const svm_node*svm_q,const svm_model*model,int dim,double sum_alpha_pos,double sum_alpha_neg,double*center_pos,double*center_neg,double& radius_pos,double& radius_neg,SVM_stat& stat)
{
	double lb_neg;
	double ub_pos;
	double L_pos;
	double U_neg;

	double lb_pos;
	double ub_neg;
	double L_neg;
	double U_pos;

	double temp_obt_dist_pos;
	double temp_obt_dist_neg;

	double LB;
	double UB;

	if(sum_alpha_pos<epsilon)
		ub_pos=0;
	else
		ub_pos=u_tri(q,center_pos,dim,radius_pos,temp_obt_dist_pos); //G_{+}

	if(sum_alpha_neg<epsilon)
		lb_neg=0;
	else
		lb_neg=ell_tri(q,center_neg,dim,radius_neg,temp_obt_dist_neg); //G_{-}

	L_pos=sum_alpha_pos*exp(-stat.gammaValue*ub_pos*ub_pos); //LB_{G_{+}}
	U_neg=sum_alpha_neg*exp(-stat.gammaValue*lb_neg*lb_neg); //UB_{G_{-}}

	//G_{+}
	if(sum_alpha_pos<epsilon)
		lb_pos=0;
	else
	{
		if(temp_obt_dist_pos<radius_pos)
			lb_pos=0;
		else
			lb_pos=temp_obt_dist_pos-radius_pos;
	}
	
	//G_{-}
	if(sum_alpha_neg<epsilon)
		ub_neg=0;
	else
		ub_neg=temp_obt_dist_neg+radius_neg;

	L_neg=sum_alpha_neg*exp(-stat.gammaValue*ub_neg*ub_neg); //LB_{G_{-}}
	U_pos=sum_alpha_pos*exp(-stat.gammaValue*lb_pos*lb_pos); //UB_{G_{+}}

	LB=L_pos-U_neg;
	UB=U_pos-L_neg;

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

	prediction(model,svm_q,stat);
}