#ifndef INIT_KC_H
#define INIT_KC_H

#include "Library.h"

//non-change constant
const double inf=9999999999999;
const double epsilon=0.000000001;

#define COMPONENT_CLOCK
#define PRUNE_STATS
//#define EXACT_VALUE_STATS

struct SVM_stat
{
	//Online for query
	double qSquareNorm;

	//Stat for the output of SVM training parameters
	double gammaValue;
	int total_sv;
	double rho;

	//Stat for the output of SVM testing
	int pruneCount;
	vector<int> class_resultVector;
	vector<int> incr_Level_Counter; //Used for Incremental Classification

	#ifdef PRUNE_STATS
		double numExact;
		double numBound;
		//The time for bound computation and exact computation
		double boundTime;
		double exactTime;
	#endif

	#ifdef COMPONENT_CLOCK
		double EC_Time;
		double BC_Time;
	#endif

	#ifdef EXACT_VALUE_STATS
		vector<double> exactValueVector;
	#endif
};

void initArray(double**& featureArray,double*& outputArray,int n,int dim);
void extractModelStatistics(fstream& file,SVM_stat& stat);
void extract_FeatureVector(char*fileName,int& n,int dim,double**& featureArray,double*& outputArray,bool isModelFile,SVM_stat& stat);
void outputResultFile(char*classResultName,SVM_stat& stat);
void outputCheck(double** featureArray,double* outputArray,int n,int dim);
void saveMatrix_toFile(char*fileName,double**matrix,int numOfElements,int dim);

#ifdef EXACT_VALUE_STATS
void output_exact_ResultFile(char*classResultName,vector<double>& exactValueVector);
#endif

#endif