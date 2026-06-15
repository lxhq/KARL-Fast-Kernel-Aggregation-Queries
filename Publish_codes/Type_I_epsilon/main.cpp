#include "GBF_KAQ.h"

int main(int argc,char**argv)
{
	char*querysetFileName=argv[1];
	char*datasetFileName=argv[2];
	char*resultFileName=argv[3];
	int method=atoi(argv[4]);
	int leafCapacity=atoi(argv[5]);
	double rel_error=atof(argv[6]);
	double b=atof(argv[7]);
	int internalCapacity=20; //only used for m-tree
	char*bulkLoad_TreeName=(char*)"";

	if(method>=4 && method<=7)
		bulkLoad_TreeName=argv[8];

	int dim;
	int n_q;
	double**queryMatrix;
	double**dataMatrix;

	KDE_stat stat={};
	stat.rel_error=rel_error;

	extract_FeatureVector(querysetFileName,n_q,dim,queryMatrix,stat);
	extract_FeatureVector(datasetFileName,stat.n,dim,dataMatrix,stat);

	preprocess_Data(dataMatrix,stat.n,dim,false,b);
	preprocess_Data(queryMatrix,n_q,dim,true,b);

	KAQ_Algorithm(queryMatrix,dataMatrix,n_q,dim,leafCapacity,internalCapacity,method,stat,bulkLoad_TreeName);

	outputResultFile(resultFileName,stat);
}
