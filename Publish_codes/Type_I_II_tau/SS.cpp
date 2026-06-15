#include "SS.h"

void SS_iter(double*q,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	double incr_Value=0;
	double sq_Value;

	for(int i=0;i<stat.total_sv;i++)
	{
		sq_Value=0;
		for(int d=0;d<dim;d++)
			sq_Value+=(q[d]-dataMatrix[i][d])*(q[d]-dataMatrix[i][d]);

		incr_Value=incr_Value+alphaArray[i]*exp(-stat.gammaValue*sq_Value);
	}

	#ifdef EXACT_VALUE_STATS
		stat.exactValueVector.push_back(incr_Value);
	#endif

	if(incr_Value>stat.rho)
		stat.class_resultVector.push_back(1);
	else
		stat.class_resultVector.push_back(-1);
}

void SS_sparse_iter(double*q,double*alphaArray,int dim,SVM_stat& stat,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,double q_sqNorm,double*p_sqNormArray,double**lu_Table)
{
	//initialization of sq_EuclidVector
	double incr_Value=0;
	double exp_LB,exp_UB;
	double incr_LB=0;
	double incr_UB=0;
	for(int i=0;i<stat.total_sv;i++)
		sq_EuclidVector[i]=q_sqNorm+p_sqNormArray[i];

	for(int d=0;d<dim;d++)
		for(int i=0;i<(int)inv_Index[d].size();i++)
			sq_EuclidVector[inv_Index[d][i].id]-=2*q[d]*inv_Index[d][i].value;

	for(int i=0;i<stat.total_sv;i++)
	{
		//incr_Value+=alphaArray[i]*exp(-stat.gammaValue*sq_EuclidVector[i]);
		lookUp_sim_Bounds(stat.gammaValue*sq_EuclidVector[i],exp_LB,exp_UB,lu_Table);
		incr_LB+=alphaArray[i]*exp_LB;
		incr_UB+=alphaArray[i]*exp_UB;
	}

	if(incr_LB>stat.rho)
		stat.class_resultVector.push_back(1);
	if(incr_UB<stat.rho)
		stat.class_resultVector.push_back(-1);
}


//This function is only used for testing
//*********************************************************************//
void convert_to_1d(double**dataMatrix,double*& dataMatrix_1d,int dim,SVM_stat& stat)
{
	dataMatrix_1d=new double[stat.total_sv*dim];
	for(int i=0;i<stat.total_sv;i++)
		for(int d=0;d<dim;d++)
			dataMatrix_1d[i*dim+d]=dataMatrix[i][d];
}
void SS_iter_1d(double*q,double*dataMatrix_1d,double*alphaArray,int dim,SVM_stat& stat)
{
	double incr_Value=0;
	double sq_Value;
	for(int i=0;i<stat.total_sv;i++)
	{
		sq_Value=0;
		for(int d=0;d<dim;d++)
			sq_Value+=(q[d]-dataMatrix_1d[i*dim+d])*(q[d]-dataMatrix_1d[i*dim+d]);

		incr_Value=incr_Value+alphaArray[i]*exp(-stat.gammaValue*sq_Value);
	}
	if(incr_Value>stat.rho)
		stat.class_resultVector.push_back(1);
	else
		stat.class_resultVector.push_back(-1);
}
//*********************************************************************//