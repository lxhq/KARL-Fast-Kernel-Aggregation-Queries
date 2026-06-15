#include "inverted_Index.h"

void build_Inverted_Index(double**dataMatrix,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,int dim,SVM_stat& stat,double*& p_sqNormArray)
{
	inv_Index=new vector<inverted_entry>[dim];
	inverted_entry tempEntry;
	p_sqNormArray=new double[stat.total_sv];
	for(int i=0;i<stat.total_sv;i++)
		p_sqNormArray[i]=0;

	for(int i=0;i<stat.total_sv;i++)
	{
		for(int d=0;d<dim;d++)
		{
			if(fabs(dataMatrix[i][d])>epsilon)
			{
				tempEntry.id=i;
				tempEntry.value=dataMatrix[i][d];
				inv_Index[d].push_back(tempEntry);
			}

			p_sqNormArray[i]+=dataMatrix[i][d]*dataMatrix[i][d];
		}
		sq_EuclidVector.push_back(0);
	}
}