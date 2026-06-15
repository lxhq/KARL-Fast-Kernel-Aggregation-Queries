#include "partition_quality_stat.h"

int main(int argc,char**argv)
{
	char*datasetFileName=argv[1];
	int method=atoi(argv[2]);
	int leafCapacity=atoi(argv[3]);
	int sampleSize=atoi(argv[4]);
	char*sampleMatrixFileName=argv[5];
	double rho=atof(argv[6]);
	double b=atof(argv[7]);

	char*bulkLoad_TreeName=(char*)"";
	int internalCapacity=20; //only used for m-tree

	int dim;
	int n;
	double**dataMatrix;
	double*alphaArray;

	SVM_stat stat;
	stat.pruneCount=0;

	if(method>=27 && method<=29) //bulk-loading methods
		bulkLoad_TreeName=argv[6];

	extract_FeatureVector_KDC(datasetFileName,n,dim,dataMatrix);
	preprocess_Data(dataMatrix,n,dim,false,b);

	update_stat(n,rho,stat);
	update_constantArray(n,alphaArray);

	offline_KC(sampleMatrixFileName,dataMatrix,alphaArray,sampleSize,dim,leafCapacity,internalCapacity,method,stat,bulkLoad_TreeName);
}