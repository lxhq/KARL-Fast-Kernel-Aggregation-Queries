//#include "partition_quality_vol.h"
#include "partition_quality_stat.h"

int main(int argc,char**argv)
{
	int dim=atoi(argv[1]);
	char*datasetFileName=argv[2];
	int method=atoi(argv[3]);
	int leafCapacity=atoi(argv[4]);
	int sampleSize=atoi(argv[5]);
	char*sampleMatrixFileName=argv[6];
	char*bulkLoad_TreeName=(char*)"";
	int internalCapacity=20; //only used for m-tree

	/*int dim=9;
	char*datasetFileName=(char*)"shuttle_model_6";
	int method=34;
	int leafCapacity=80;
	int sampleSize=100;
	char*sampleMatrixFileName=(char*)"shuttle_dataSample_100";
	char*bulkLoad_TreeName=(char*)"";*/

	int n;
	double**dataMatrix;
	double*alphaArray;

	SVM_stat stat;
	stat.pruneCount=0;

	if(method>=27 && method<=29) //bulk-loading methods
		bulkLoad_TreeName=argv[7];

	extract_FeatureVector(datasetFileName,n,dim,dataMatrix,alphaArray,true,stat);	

	offline_KC(sampleMatrixFileName,dataMatrix,alphaArray,sampleSize,dim,leafCapacity,internalCapacity,method,stat,bulkLoad_TreeName);

	//quality_Check(dataMatrix,alphaArray,dim,leafCapacity,internalCapacity,method,stat,bulkLoad_TreeName);
}