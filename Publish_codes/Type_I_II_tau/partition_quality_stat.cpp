#include "partition_quality_stat.h"

void sample_Dataset(char*sampleMatrixFileName,double**dataMatrix,double**& sampleMatrix,int dataSize, int sampleSize,int dim)
{
	fstream sampleFile;
	fstream newSampleFile;
	int id;
	sampleMatrix=new double*[sampleSize];
	for(int s=0;s<sampleSize;s++)
		sampleMatrix[s]=new double[dim];

	sampleFile.open(sampleMatrixFileName);
	if(sampleFile.is_open()==true)
	{
		for(int s=0;s<sampleSize;s++)
			for(int d=0;d<dim;d++)
				sampleFile>>sampleMatrix[s][d];

		sampleFile.close();
		return;
	}

	cout<<"No Sample Matrix!"<<endl;
	srand(time(NULL));

	for(int s=0;s<sampleSize;s++)
	{
		id=rand()%dataSize;
		for(int d=0;d<dim;d++)
			sampleMatrix[s][d]=dataMatrix[id][d];
	}

	newSampleFile.open(sampleMatrixFileName,ios::in | ios::out | ios::trunc);
	if(newSampleFile.is_open()==false)
	{
		cout<<"Error in saving sampleFile!"<<endl;
		exit(0);
	}

	for(int s=0;s<sampleSize;s++)
	{
		for(int d=0;d<dim;d++)
			newSampleFile<<sampleMatrix[s][d]<<" ";
		newSampleFile<<endl;
	}

	newSampleFile.close();
	cout<<"Finish saving the sample file!"<<endl;
}

void offline_KC(char*sampleMatrixFileName,double**dataMatrix,double*alphaArray,int sampleSize,int dim,int leafCapacity,int internalCapacity,int method,SVM_stat& stat,char*bulkLoad_TreeName)
{
	double**sampleMatrix;
	sample_Dataset(sampleMatrixFileName,dataMatrix,sampleMatrix,stat.total_sv,sampleSize,dim);

	kdTree kd_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);
	binaryTree binary_Tree(dim,dataMatrix,alphaArray,leafCapacity,stat);
	mTree m_Tree(dim,dataMatrix,alphaArray,internalCapacity,leafCapacity,stat);

	//Used in LibSVM linear scan method

	svm_node**svm_qMatrix;
	convert_2D_qMatrix_to_SVMFormat(sampleMatrix,svm_qMatrix,dim,sampleSize);

	if(method==1)//tKDC kd-tree + LB_MBR and UB_MBR
		kd_Tree.rootNode=new kdNode();
	if(method==4)//kd-tree + LinearApprox
		kd_Tree.rootNode=new kdLinearAugNode();
	if(method==5)//binary-tree + LinearApprox + GPA
		binary_Tree.rootNode=new binaryNode();
	if(method==27)//m-tree + triangle inequality
		m_Tree.rootNode=new mNode();
	if(method==28)//m-tree + LB_MBR and UB_MBR
		m_Tree.rootNode=new mAugNode();
	if(method==29)//m-tree + LinearApprox
		m_Tree.rootNode=new mLinearAugNode();
	if(method==24) //binary-tree + LinearApprox + GPA (norm partition)
		binary_Tree.rootNode=new normNode();
	if(method==25)
		binary_Tree.rootNode=new farNode();
	if(method==26)
		binary_Tree.rootNode=new ipNode();
	if(method==34)
		binary_Tree.rootNode=new ballNode();
	if(method==35)
		binary_Tree.rootNode=new ballNode_SOTA();

	//build tree Preprocessing + create augment tree
	if(method>=1 && method<=4)
	{
		kd_Tree.build_kdTree(stat);
		kd_Tree.updateAugment((kdNode*)kd_Tree.rootNode);
	}

	if(method==5 || method==24 || method==25 || method==26 || method==34 || method==35) //binary tree
		binary_Tree.build_BinaryTree();
		
	if(method==27 || method==28 || method==29)
		m_Tree.load_Tree(bulkLoad_TreeName);

	#ifdef PRUNE_STATS
		stat.numBound=0;
		stat.numExact=0;
		stat.boundTime=0;
		stat.exactTime=0;
		binary_Tree.stat.numBound=0;
		binary_Tree.stat.numExact=0;
		binary_Tree.stat.boundTime=0;
		binary_Tree.stat.exactTime=0;
	#endif


	#ifdef C_PLUSPLUS11_CLOCK
		auto start_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_s=clock();
	#endif

	for(int q=0;q<sampleSize;q++)
	{
		//Online One-Time processing
		switch (method)
		{
			case 1:
				GBF_iter(sampleMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 4:
				stat.qSquareNorm=computeSqNorm(sampleMatrix[q],dim);
				GBF_iter(sampleMatrix[q],svm_qMatrix[q],kd_Tree,dim,stat);
				break;
			case 5:
			case 24:
			case 25:
			case 26:
			case 34:
			case 35:
				stat.qSquareNorm=computeSqNorm(sampleMatrix[q],dim);
				GBF_iter(sampleMatrix[q],svm_qMatrix[q],binary_Tree,dim,stat);
				break;
			case 27:
			case 28:
				GBF_iter(sampleMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
			case 29:
				stat.qSquareNorm=computeSqNorm(sampleMatrix[q],dim);
				GBF_iter(sampleMatrix[q],svm_qMatrix[q],m_Tree,dim,stat);
				break;
		}
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_s=chrono::high_resolution_clock::now();
		double online_Time=(chrono::duration_cast<chrono::nanoseconds>(end_s-start_s).count())/1000000000.0;
	#else
		clock_t end_s=clock();
		double online_Time=((double)(end_s-start_s))/CLOCKS_PER_SEC;
	#endif

	cout<<"method "<<method<<": "<<sampleSize/online_Time<<"samples/sec"<<endl;
}
