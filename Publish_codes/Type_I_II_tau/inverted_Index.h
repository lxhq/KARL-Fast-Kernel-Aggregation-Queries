#ifndef INVERTED_INDEX_H
#define INVERTED_INDEX_H

#include "init_KC.h"
struct inverted_entry
{
	int id;
	double value;
};

void build_Inverted_Index(double**dataMatrix,vector<inverted_entry>*& inv_Index,vector<double>& sq_EuclidVector,int dim,SVM_stat& stat,double*& p_sqNormArray);

#endif