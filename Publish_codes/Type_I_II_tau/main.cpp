#include "init_KC.h"
#include "GBF_KC.h"

int main(int argc,char**argv)
{
	int dim=atoi(argv[1]);
	int qNum=atoi(argv[2]);
	char*querysetFileName=argv[3];
	char*datasetFileName=argv[4];
	char*classResultName=argv[5];
	int method=atoi(argv[6]);
	int leafCapacity=atoi(argv[7]);
	int internalCapacity=20; //only used for m-tree

	char*oracleFileName;
	char*pivotFileName;
	char*bulkLoad_TreeName;

	int n;
	double**queryMatrix;
	double*queryOutputArray;
	double**dataMatrix;
	double*alphaArray;	

	//Used by libSVM method
	svm_model*model;

	SVM_stat stat;

	stat.pruneCount=0;

	extract_FeatureVector(querysetFileName,qNum,dim,queryMatrix,queryOutputArray,false,stat);
	extract_FeatureVector(datasetFileName,n,dim,dataMatrix,alphaArray,true,stat);

	//default chunkSize
	int chunkSize=(int)ceil(stat.total_sv/10.0);

	if(method==19) //Statistically obtain the chunkSize
		chunkSize=loadChunkSize(argv[8]);
	if(method==100 || method==101 || method==102) //Oracle 
		oracleFileName=argv[8];
	if(method==20 || method==22) //Select Pivots and Reorder
		pivotFileName=argv[8];
	if(method>=27 && method<=29) //bulk-loading methods
		bulkLoad_TreeName=argv[8];
	
	KC_Algorithm(queryMatrix,dataMatrix,alphaArray,qNum,dim,leafCapacity,internalCapacity,method,stat,model,chunkSize,oracleFileName,pivotFileName,bulkLoad_TreeName);

	#ifdef EXACT_VALUE_STATS
		output_exact_ResultFile(classResultName,stat.exactValueVector);
		return 0;
	#endif

	outputResultFile(classResultName,stat);
}

//if(method==13 || method==16) //input for testing the mnist dataset
//	chunkSize=atoi(argv[8]);
