#include "Reorder.h"

//sort in increasing order
bool compare_order(const orderEntry& e1,const orderEntry& e2)
{
	return (e1.dist<e2.dist);
}

bool compare_order_norm(const orderEntry& e1,const orderEntry& e2)
{
	return (e1.norm<e2.norm);
}

void concentric_Reorder(double**P,vector<orderEntry>& orderList,double*alphaArray,int dim,SVM_stat& stat)
{
	double*centroid=new double[dim];
	double sum_alpha=0;
	orderEntry tempEntry;

	for(int d=0;d<dim;d++)
		centroid[d]=0;

	for(int i=0;i<stat.total_sv;i++)
	{
		sum_alpha+=alphaArray[i];
		for(int d=0;d<dim;d++)
			centroid[d]+=alphaArray[i]*P[i][d];
	}

	for(int d=0;d<dim;d++)
		centroid[d]=centroid[d]/sum_alpha;

	//Compute the distance from centroid to all objects
	for(int i=0;i<stat.total_sv;i++)
	{
		tempEntry.id=i;
		tempEntry.dist=euclid_dist(centroid,P[i],dim);
		orderList.push_back(tempEntry);
	}

	sort(orderList.begin(),orderList.end(),compare_order);
}

void norm_Reorder(double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat)
{
	double sq_norm;
	orderEntry tempEntry;
	for(int i=0;i<stat.total_sv;i++)
	{
		sq_norm=0;
		for(int d=0;d<dim;d++)
			sq_norm+=P[i][d]*P[i][d];

		tempEntry.id=i;
		tempEntry.norm=sq_norm;

		orderList.push_back(tempEntry);
	}

	sort(orderList.begin(),orderList.end(),compare_order_norm);
}

void reference_Reorder(double*q,double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat)
{

}

void oracle_Reorder(double*q,double**P,vector<orderEntry>& orderList,int dim,SVM_stat& stat)
{
	orderEntry tempEntry;
	for(int i=0;i<stat.total_sv;i++)
	{
		tempEntry.id=i;
		tempEntry.dist=euclid_dist(q,P[i],dim);
		orderList.push_back(tempEntry);
	}

	sort(orderList.begin(),orderList.end(),compare_order);
}

void reordering(double**P,double*alphaArray,vector<orderEntry>& orderList,SVM_stat& stat,int dim)
{
	double**PP=new double*[stat.total_sv];
	double*reorder_alphaArray=new double[stat.total_sv];

	for(int i=0;i<stat.total_sv;i++)
		PP[i]=new double[dim];

	for(int i=0;i<stat.total_sv;i++)
	{
		reorder_alphaArray[i]=alphaArray[orderList[i].id];
		for(int d=0;d<dim;d++)
			PP[i][d]=P[orderList[i].id][d];
	}

	//Copy back to P and alphaArray
	for(int i=0;i<stat.total_sv;i++)
	{
		alphaArray[i]=reorder_alphaArray[i];
		for(int d=0;d<dim;d++)
			P[i][d]=PP[i][d];
	}
	
	delete[] reorder_alphaArray;
	for(int i=0;i<stat.total_sv;i++)
		delete[] PP[i];
	delete[] PP;
}

void reordering_mless(double**P,double**PP,double*alphaArray,double*reorder_alphaArray,vector<orderEntry>& orderList,SVM_stat& stat,int dim)
{
	for(int i=0;i<stat.total_sv;i++)
	{
		reorder_alphaArray[i]=alphaArray[orderList[i].id];
		for(int d=0;d<dim;d++)
			PP[i][d]=P[orderList[i].id][d];
	}

	//Copy back to P and alphaArray
	for(int i=0;i<stat.total_sv;i++)
	{
		alphaArray[i]=reorder_alphaArray[i];
		for(int d=0;d<dim;d++)
			P[i][d]=PP[i][d];
	}
}

void init_reorder_mless(double**& PP,double*& reorder_alphaArray,SVM_stat& stat,int dim)
{
	PP=new double*[stat.total_sv];
	reorder_alphaArray=new double[stat.total_sv];

	for(int i=0;i<stat.total_sv;i++)
		PP[i]=new double[dim];
}

void delete_reorder_mless(double**PP,double*reorder_alphaArray,SVM_stat& stat)
{
	for(int i=0;i<stat.total_sv;i++)
		delete[] PP[i];

	delete[] PP;
	delete[] reorder_alphaArray;
}