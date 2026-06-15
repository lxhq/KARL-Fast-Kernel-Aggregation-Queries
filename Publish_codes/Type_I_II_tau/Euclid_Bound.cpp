#include "Euclid_Bound.h"

double ell_MBR(double*q,double**boundary,int dim)
{
	double ell=0;
	double temp;
	for(int d=0;d<dim;d++)
	{
		if(q[d]<boundary[d][0])
			temp=boundary[d][0]-q[d];
		else
		{
			if(q[d]>boundary[d][1])
				temp=q[d]-boundary[d][1];
			else
				temp=0;
		}
		
		ell+=temp*temp;
	}
	ell=sqrt(ell);

	return ell;
}

double u_MBR(double*q,double**boundary,int dim)
{
	double u=0;
	double temp;
	for(int d=0;d<dim;d++)
	{
		temp=max(fabs(q[d]-boundary[d][0]),fabs(q[d]-boundary[d][1]));
		u+=temp*temp;
	}
	u=sqrt(u);

	return u;
}

double oneD_nearest(double v,double*Array,int left,int right)
{
	int middle;
	double diff;
	if(left==right) //only one element
		return Array[left];
	if(right-left==1) //two elements
		return min(fabs(v-Array[left]),fabs(Array[right]-v));

	middle=(int)floor((double)(left+right)/2);

	diff=abs(v-Array[middle])-abs(v-Array[middle+1]);

	if(diff>epsilon)
		return oneD_nearest(v,Array,middle+1,right);
	else
	{
		if(diff>=-epsilon && diff<=epsilon)
			return abs(v-Array[middle]);
		else //diff<-epsilon
			return oneD_nearest(v,Array,left,middle);
	}	
}

double ell_MIN(double*q,double**boundary,double**d_SortObj,int lastIndex,int dim)
{
	double ell=0;
	double temp;
	for(int d=0;d<dim;d++)
	{
		if(q[d]>boundary[d][0] && q[d]<boundary[d][1])
			temp=oneD_nearest(q[d],d_SortObj[d],0,lastIndex);

		if(q[d]<=boundary[d][0])
			temp=boundary[d][0]-q[d];
		if(q[d]>=boundary[d][1])
			temp=q[d]-boundary[d][1];
		
		ell+=temp*temp;
	}
	ell=sqrt(ell);

	return ell;
}

double u_tri(double*q,double*center,int dim,double radius,double& obt_dist)
{
	double u;
	obt_dist=0;

	for(int d=0;d<dim;d++)
		obt_dist+=(q[d]-center[d])*(q[d]-center[d]);

	obt_dist=sqrt(obt_dist);

	u=obt_dist+radius;

	return u;

}

double euclid_dist(double*q,double*p,int dim)
{
	double dist=0;
	for(int d=0;d<dim;d++)
		dist+=(q[d]-p[d])*(q[d]-p[d]);
	return sqrt(dist);
}