#ifndef EUCLID_BOUND_H
#define EUCLID_BOUND_H

#include "init_KC.h"

//No extern--> cannot compile
double ell_MBR(double*q,double**boundary,int dim);
double u_MBR(double*q,double**boundary,int dim);
double ell_tri(double*q,double*center,int dim,double radius,double& obt_dist);
double u_tri(double*q,double*center,int dim,double radius,double& obt_dist);
double euclid_dist(double*q,double*p,int dim);

#endif