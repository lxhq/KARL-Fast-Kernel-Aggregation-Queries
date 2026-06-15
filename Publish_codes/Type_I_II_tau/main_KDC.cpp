#include "init_KC.h"
#include "GBF_KC.h"

int main(int argc,char**argv)
{
	char*querysetFileName=argv[1];
	char*datasetFileName=argv[2];
	char*classResultName=argv[3];
	int method=atoi(argv[4]);
	int leafCapacity=atoi(argv[5]);
	double rho=atof(argv[6]);
	double b=atof(argv[7]);
	int internalCapacity=20; //only used for m-tree

	char*bulkLoad_TreeName=(char*)"";

	int qNum;
	int pNum;
	int dim;
	double**queryMatrix;
	double**dataMatrix;
	double*alphaArray;

	char*notUsed=(char*)"";

	//Used by libSVM method
	svm_model*model;
	SVM_stat stat;
	stat.pruneCount=0;

	extract_FeatureVector_KDC(datasetFileName,pNum,dim,dataMatrix);
	extract_FeatureVector_KDC(querysetFileName,qNum,dim,queryMatrix);

	preprocess_Data(dataMatrix,pNum,dim,false,b);
	preprocess_Data(queryMatrix,qNum,dim,true,b);

	update_stat(pNum,rho,stat);
	update_constantArray(pNum,alphaArray);

	if(method>=27 && method<=29) //bulk-loading methods
		bulkLoad_TreeName=argv[8];
	
	KC_Algorithm(queryMatrix,dataMatrix,alphaArray,qNum,dim,leafCapacity,internalCapacity,method,stat,model,0,notUsed,notUsed,bulkLoad_TreeName);

	#ifdef EXACT_VALUE_STATS
		output_exact_ResultFile(classResultName,stat.exactValueVector);
		return 0;
	#endif

	outputResultFile(classResultName,stat);
}
