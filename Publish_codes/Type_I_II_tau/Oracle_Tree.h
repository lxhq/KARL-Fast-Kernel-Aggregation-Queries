#ifndef ORACLE_TREE_H
#define ORACLE_TREE_H

#include "GBF_KC.h"
#include "binaryTree.h"
#include "Reorder.h"

double runOracleTree(double*q,svm_node*svm_q,binaryTree& oracle_Tree);
void runOracleTree_mless(double*q,svm_node*svm_q,double**& PP,double*& reorder_alphaArray,binaryTree& oracle_Tree,double& online_Time);

#endif