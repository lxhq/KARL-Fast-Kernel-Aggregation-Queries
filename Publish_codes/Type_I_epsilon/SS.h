#ifndef SS_H
#define SS_H

#include "init_KAQ.h"
#include "Validation.h"
#include "Euclid_Bound.h"

void SS_iter(double*q,double**dataMatrix,int dim,KDE_stat& stat);

//method=8
void pre_Compute_sequential_MBR(double**dataMatrix,int dim,double**& boundary,KDE_stat& stat);
void SS_MBR(double*q,double**dataMatrix,int dim,double**boundary,KDE_stat& stat);

//method=9
void pre_Compute_sequential_Delta(double**dataMatrix,int dim,double*& center,double& radius,KDE_stat& stat);
void SS_Delta(double*q,double**dataMatrix,int dim,double*center,double radius,KDE_stat& stat);

//method=10
void pre_Compute_sequential(double**dataMatrix,int dim,double**& boundary,double*& a_G,double& S_G,KDE_stat& stat);
void SS_linear(double*q,double**dataMatrix,int dim,double**& boundary,double*a_G,double S_G,KDE_stat& stat);
#endif