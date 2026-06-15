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

#ifdef KARL_PROFILE
	long long profile_queries;
	long long profile_nodes_processed;
	long long profile_leaf_nodes;
	long long profile_exact_points;
	long long profile_bound_calls;
	long long profile_validate_checks;
	long long profile_heap_pushes;
	long long profile_heap_pops;
	long long profile_validate_success_queries;
	long long profile_exact_finish_queries;
	double profile_time_validate_sec;
	double profile_time_bound_sec;
	double profile_time_leaf_eval_sec;
	vector<long long> profile_nodes_per_query;
	vector<long long> profile_leaf_nodes_per_query;
	vector<long long> profile_exact_points_per_query;
#endif
};

void initArray(double**& featureArray,int n,int dim);
void extract_FeatureVector(char*fileName,int& n,int& dim,double**& featureArray,KDE_stat& stat);
void preprocess_Data(double**featureArray,int n,int dim,bool isQuery,double b);
void outputResultFile(char*resultFileName,KDE_stat& stat);

#endif
