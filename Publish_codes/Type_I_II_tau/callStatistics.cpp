#include "init_KC.h"
#include "GBF_KC.h"
#include "statistics.h"

int main(int argc,char**argv)
{
	int dim=atoi(argv[1]);
	int qNum=atoi(argv[2]);
	char*querysetFileName=argv[3];
	char*datasetFileName=argv[4];
	int method=atoi(argv[5]);
	int leafCapacity=atoi(argv[6]);
	int module=atoi(argv[7]);
	int model_id=atoi(argv[8]);
	double rho;
	double b;
	

	//debug
	/*int dim=300;
	int qNum=10000;
	char*querysetFileName=(char*)"w8a_sample.t";
	char*datasetFileName=(char*)"w8a_model_validate";
	int method=1;
	int leafCapacity=20;
	int module=1;
	int model_id=3;
	double rho;
	double b;*/

	int internalCapacity=20; //only used for m-tree
	char*bulkLoad_TreeName=(char*)"";

	int n;
	double**queryMatrix;
	double*queryOutputArray;
	double**dataMatrix;
	double*alphaArray;
	bool is_plus; //used in model 3

	SVM_stat stat;

	stat.pruneCount=0;

	if(model_id==1) //KDC
	{
		rho=atof(argv[9]);
		b=atof(argv[10]);
		
		if(method>=27 && method<=29) //bulk-loading methods
			bulkLoad_TreeName=argv[11];
		extract_FeatureVector_KDC(datasetFileName,n,dim,dataMatrix);
		extract_FeatureVector_KDC(querysetFileName,qNum,dim,queryMatrix);

		preprocess_Data(dataMatrix,n,dim,false,b);
		preprocess_Data(queryMatrix,qNum,dim,true,b);

		update_stat(n,rho,stat);
		update_constantArray(n,alphaArray);
	}
	if(model_id==2) //OCSVM
	{
		extract_FeatureVector(querysetFileName,qNum,dim,queryMatrix,queryOutputArray,false,stat);
		extract_FeatureVector(datasetFileName,n,dim,dataMatrix,alphaArray,true,stat);
		if(method>=27 && method<=29) //bulk-loading methods
			bulkLoad_TreeName=argv[9];
	}
	if(model_id==3) //SVM
	{
		is_plus=(bool)atoi(argv[9]);
		extract_FeatureVector_twoClass_SVM_stat(querysetFileName,qNum,dim,queryMatrix,queryOutputArray,false,stat,is_plus);
		extract_FeatureVector_twoClass_SVM_stat(datasetFileName,n,dim,dataMatrix,alphaArray,true,stat,is_plus);
		if(method>=27 && method<=29) //bulk-loading methods
			bulkLoad_TreeName=argv[10];
	}
	
	//default chunkSize
	int chunkSize=(int)ceil(stat.total_sv/10.0);

	cout<<"Model: "<<model_id<<endl;
	cout<<"Dataset: "<<datasetFileName<<endl;
	cout<<"Module:"<<module<<endl;

	if(model_id==3)
	{
		if(is_plus==true)
			cout<<"Support Vectors (Positive)"<<endl;
		else
			cout<<"Support Vectors (negative)"<<endl;
	}
	statModule_Selection(queryMatrix,dataMatrix,alphaArray,qNum,dim,leafCapacity,internalCapacity,stat,method,module,bulkLoad_TreeName);
}