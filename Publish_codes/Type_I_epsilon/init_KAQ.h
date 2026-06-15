#ifndef INIT_KAQ_H
#define INIT_KAQ_H

#include "Library.h"

//non-change constant
const double inf=9999999999999;
const double epsilon=0.000000001;

struct KDE_stat
{
	//Online for query
	double qSquareNorm;

	//Stat for KDE
	int n;
	double rel_error;

	//Output result value vectors
	vector<double> resultValueVector;
};

void initArray(double**& featureArray,int n,int dim);
void extract_FeatureVector(char*fileName,int& n,int& dim,double**& featureArray,KDE_stat& stat);
void preprocess_Data(double**featureArray,int n,int dim,bool isQuery,double b);
void outputResultFile(char*resultFileName,KDE_stat& stat);

#endif