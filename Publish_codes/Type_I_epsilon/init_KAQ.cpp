#include "init_KAQ.h"

void initArray(double**& featureArray,int n,int dim)
{
	featureArray=new double*[n];
	for(int i=0;i<n;i++)
		featureArray[i]=new double[dim];

	//initalization of featureArray
	for(int i=0;i<n;i++)
		for(int d=0;d<dim;d++)
			featureArray[i][d]=0;
}

void extract_FeatureVector(char*fileName,int& n,int& dim,double**& featureArray,KDE_stat& stat)
{
	fstream file;
	file.open(fileName);
	if(file.is_open()==false)
	{
		cout<<"Error in opening the file!"<<endl;
		return;
	}

	file>>n;
	file>>dim;

	initArray(featureArray,n,dim);

	for(int i=0;i<n;i++)
		for(int d=0;d<dim;d++)
			file>>featureArray[i][d];
}

double obtain_std(double**featureArray,int n,int d)
{
	double sq_sum=0;
	double sum=0;
	double mean;
	double std;

	for(int i=0;i<n;i++)
	{
		sq_sum+=(featureArray[i][d])*(featureArray[i][d]);
		sum+=featureArray[i][d];
	}
	mean=sum/n;
	std=sqrt((sq_sum-n*mean*mean)/n);

	return std;
}

void preprocess_Data(double**featureArray,int n,int dim,bool isQuery,double b)
{
	//double b=10;
	//double constant=pow((double)n,1.0/((double)(dim+4)))/2.0;
	double constant=pow((double)n,1.0/((double)(dim+4)))/(2.0*b);
	static double*coeff_vec;
	double std;

	if(isQuery==false)
	{
		coeff_vec=new double[dim];

		//obtain the coeff_vec
		for(int d=0;d<dim;d++)
		{
			std=obtain_std(featureArray,n,d);
			coeff_vec[d]=sqrt(constant/std);
		}
	}

	for(int i=0;i<n;i++)
		for(int d=0;d<dim;d++)
			featureArray[i][d]*=coeff_vec[d];
}

void outputResultFile(char*resultFileName,KDE_stat& stat)
{
	fstream resultFile;
	resultFile.open(resultFileName,ios::in|ios::out|ios::trunc);

	if(resultFile.is_open()==false)
	{
		cout<<"Cannot Open Result File!"<<endl;
		return;
	}

	for(int r=0;r<(int)stat.resultValueVector.size();r++)
		resultFile<<stat.resultValueVector[r]<<endl;

	resultFile.close();
}