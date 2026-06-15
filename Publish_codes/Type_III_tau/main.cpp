#include "GBF_KC_Binary.h"

int main(int argc,char**argv)
{
	int dim=atoi(argv[1]);
	int qNum=atoi(argv[2]);
	char*querysetFileName=argv[3];
	char*datasetFileName=argv[4];
	char*classResultName=argv[5];
	int method=atoi(argv[6]);
	int leafCapacity=atoi(argv[7]);
	int internalCapacity=2; //only used for m-tree

	char*bulkLoad_TreeName=(char*)"";

	if(method>=4 && method<=6)
		bulkLoad_TreeName=argv[8];

	int n;
	double**queryMatrix;
	double*queryOutputArray;
	double**dataMatrix;
	double*alphaArray;

	SVM_stat stat;

	stat.pruneCount=0;

	extract_FeatureVector(querysetFileName,qNum,dim,queryMatrix,queryOutputArray,false,stat);
	extract_FeatureVector(datasetFileName,n,dim,dataMatrix,alphaArray,true,stat);

	KC_Algorithm(queryMatrix,dataMatrix,alphaArray,qNum,dim,leafCapacity,internalCapacity,method,stat,bulkLoad_TreeName);

	#ifdef EXACT_VALUE_STATS
		output_exact_ResultFile(classResultName,stat.exactValueVector);
		return 0;
	#endif

	outputResultFile(classResultName,stat);
}
