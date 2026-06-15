#include "init_KC.h"
#include "GBF_KC.h"
#include "my_svm.h"

/*void loadData(char*datasetFileName,double**& dataMatrix,double*& alphaArray,int dim)
{
	int total_sv;
	//double alpha;
	fstream file;
	file.open(datasetFileName);
	if(file.is_open()==false)
	{
		cout<<"Cannot Open File!"<<endl;
		return;
	}
	file>>total_sv;
	file>>dim;

	dataMatrix=new double*[total_sv];
	alphaArray=new double[total_sv];
	for(int i=0;i<total_sv;i++)
		dataMatrix[i]=new double[dim];

	for(int i=0;i<total_sv;i++)
	{
		for(int d=0;d<dim;d++)
			file>>dataMatrix[i][d];
		file>>alphaArray[i];
	}
	
	file.close();
}

void loadQuery(char*queryFileName,double**& queryMatrix,int dim)
{
	int qNum;
	fstream qFile;
	qFile.open(queryFileName);
	if(qFile.is_open()==false)
	{
		cout<<"Cannot Open File!"<<endl;
		return;
	}
	qFile>>qNum;
	qFile>>dim;

	queryMatrix=new double*[qNum];
	for(int i=0;i<qNum;i++)
		queryMatrix[i]=new double[dim];

	for(int i=0;i<qNum;i++)
	{
		for(int d=0;d<dim;d++)
			qFile>>queryMatrix[i][d];
	}
	
	qFile.close();
}*/

//test tree
//int main()
//{
	/*int dim=2;
	char*queryFileName=(char*)"query_small.txt";
	char*datasetFileName=(char*)"test_TreeFile.txt";
	int method=5;
	int leafCapacity=2;
	int internalCapacity=2;

	double**queryMatrix;
	double**dataMatrix;
	double*alphaArray;
	SVM_stat stat;
	int qNum=5;
	stat.total_sv=14;

	loadQuery(queryFileName,queryMatrix,dim);
	loadData(datasetFileName,dataMatrix,alphaArray,dim);
	stat.rho=1.4;
	stat.gammaValue=1.0/stat.total_sv;

	mTree m_Tree(dim,dataMatrix,alphaArray,internalCapacity,leafCapacity,stat);
	m_Tree.rootNode=new mNode();
	((mNode*)m_Tree.rootNode)->O_r=new double[dim];

	m_Tree.build_m_tree();

 	//m_Tree.lookUp_m_tree();

	for(int q=0;q<qNum;q++)
	{
		//SS_iter(queryMatrix[q],dataMatrix,alphaArray,dim,stat);
		stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
		GBF_iter(queryMatrix[q],m_Tree,dim,stat);
	}

	for(int q=0;q<qNum;q++)
		cout<<stat.class_resultVector[q]<<" ";*/

	/*binaryTree binary_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);

	binary_Tree.rootNode=new binaryNode();

	binary_Tree.build_BinaryTree();

	binary_Tree.lookUp_BinaryTree();

	for(int q=0;q<qNum;q++)
	{
		//SS_iter(queryMatrix[q],dataMatrix,alphaArray,dim,stat);
		stat.qSquareNorm=computeSqNorm(queryMatrix[q],dim);
		GBF_iter(queryMatrix[q],binary_Tree,dim,stat);
	}

	for(int q=0;q<qNum;q++)
		cout<<stat.class_resultVector[q]<<" ";

	cout<<endl;*/

	/*mTree m_Tree(dim,dataMatrix,alphaArray,internalCapacity,leafCapacity,stat);
	m_Tree.rootNode=new mNode();

	m_Tree.build_m_tree(stat);

	m_Tree.lookUp_mTree();*/
//}

/*int main(int argc,char**argv)
{
	int dim=41;
	int qNum=1000;
	char*querysetFileName=(char*)"kddCup_q_scale.dat";
	//char*datasetFileName=(char*)"KDDCup_model";
	char*datasetFileName=(char*)"KDDCup_model_2";
	char*classResultName=(char*)"test.txt";
	int method=23;
	int leafCapacity=20;
	int internalCapacity=20; //only used for m-tree
	char*oracleFileName=(char*)"testingOracle.txt";
	char*pivotFileName=(char*)"nil";

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

	if(method==100 || method==101 || method==102)
		oracleFileName=(char*)"testing.txt";
	if(method==20 || method==22)
		pivotFileName=(char*)"KDDCup_model_2_kMean_10";

	KC_Algorithm(queryMatrix,dataMatrix,alphaArray,qNum,dim,leafCapacity,internalCapacity,method,stat,model,chunkSize,oracleFileName,pivotFileName);

	outputResultFile(classResultName,stat);
}*/