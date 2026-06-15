#ifndef EXP_LOOKUP_H
#define EXP_LOOKUP_H

#include <cmath>

const int TableSize=10000000;

//Create the lookup table for exp(-x) LB<=exp(-x)<=UB
void lookUpTable_Creation(double**& lu_Table);

//Online Computation
double lookUpLB(double x,double**lu_Table);
double lookUpUB(double x,double**lu_Table);
void lookUp_sim_Bounds(double x,double& l,double& u,double**lu_Table);

#endif