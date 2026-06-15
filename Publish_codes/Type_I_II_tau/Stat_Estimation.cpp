#include "Stat_Estimation.h"

//Offline Estimation
void stat_Estimation(char*speedUpFileName,int trainNum,double**trainMatrix,double**P,double*alphaArray,SVM_stat& stat,int dim,svm_model*& model, int PT)
{
	svm_node**svm_trainMatrix;
	group_Info*groupArray;
	vector<orderEntry> orderList;
	fstream speedUpModel;
	svm_node*x_space;

	speedUpModel.open(speedUpFileName,ios::out | ios::app);

	int chunkSize=(int)floor(stat.total_sv/(double)PT);

	convert_2D_qMatrix_to_SVMFormat(trainMatrix,svm_trainMatrix,dim,trainNum);
	init_model_inMemory(P,alphaArray,dim,stat,model,x_space);
	concentric_Reorder(P,orderList,alphaArray,dim,stat);
	reordering(P,alphaArray,orderList,stat,dim);
	create_groupArray(P,alphaArray,groupArray,chunkSize,dim,stat,true);

	for(int t=0;t<trainNum;t++)
	{
		stat.qSquareNorm=computeSqNorm(trainMatrix[t],dim);
		incremental_Classification_multiStep(trainMatrix[t],P,groupArray,alphaArray,chunkSize,dim,stat,model,svm_trainMatrix[t]);
	}

	//handle the statistics (Code here)
	double rho=(double)stat.incr_Level_Counter[PT]/trainNum;
	double rho_i;
	double denominator=(rho+(1-rho)*PT);
	double speedUp;
	for(int i=0;i<PT;i++)
	{
		rho_i=((double)stat.incr_Level_Counter[i]/trainNum);
		denominator+=((double)stat.total_sv/PT)*(i+1)*rho_i;
	}
	speedUp=(double)stat.total_sv/denominator;

	speedUpModel<<PT<<" "<<speedUp<<endl;

	speedUpModel.close();
}

void load_sample_idFile(char*sample_idFileName,double**& trainMatrix,int& tNum,double**dataMatrix,int dim)
{
	fstream sample_idFile;
	//int sampleNum;
	int id;
	sample_idFile.open(sample_idFileName);

	if(sample_idFile.is_open()==false)
	{
		cout<<"Cannot Open sample_idFile!"<<endl;
		exit(1);
	}

	sample_idFile>>tNum;

	trainMatrix=new double*[tNum];
	for(int t=0;t<tNum;t++)
		trainMatrix[t]=new double[dim];

	for(int t=0;t<tNum;t++)
	{
		sample_idFile>>id;
		for(int d=0;d<dim;d++)
			trainMatrix[t][d]=dataMatrix[id][d];
	}

	sample_idFile.close();
}

//Online Estimation
int loadChunkSize(char*speedUpFileName)
{
	int numberOfChunk;
	int chunkSize,bestChunkSize;
	double speedUpValue,highestSpeedUpValue;

	fstream speedUpFile;
	speedUpFile.open(speedUpFileName);
	if(speedUpFile.is_open()==false)
	{
		cout<<"Cannot Open speedUpFile!"<<endl;
		exit(1);
	}
	speedUpFile>>numberOfChunk;

	bestChunkSize=-1;
	highestSpeedUpValue=-1;
	for(int c=0;c<numberOfChunk;c++)
	{
		speedUpFile>>chunkSize;
		speedUpFile>>speedUpValue;

		if(speedUpValue>highestSpeedUpValue)
		{
			highestSpeedUpValue=speedUpValue;
			bestChunkSize=chunkSize;
		}
	}

	return bestChunkSize;
}