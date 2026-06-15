#include "init_KC.h"

void initArray(double**& featureArray,double*& outputArray,int n,int dim)
{
	featureArray=new double*[n];
	for(int i=0;i<n;i++)
		featureArray[i]=new double[dim];

	//initalization of featureArray
	for(int i=0;i<n;i++)
		for(int d=0;d<dim;d++)
			featureArray[i][d]=0;

	outputArray=new double[n];
}

void extractModelStatistics(fstream& file,SVM_stat& stat)
{
	string temp;
	string svm_type;
	string kernel_type;
	int nr_class;
	int*label;
	int*nr_sv;

	file>>temp;
	file>>svm_type;

	file>>temp;
	file>>kernel_type;

	file>>temp;
	file>>stat.gammaValue;

	file>>temp;
	file>>nr_class;

	//create the array
	label=new int[nr_class];
	nr_sv=new int[nr_class];

	file>>temp;
	file>>stat.total_sv;

	file>>temp;
	file>>stat.rho;

	file>>temp;file>>temp;file>>temp;
	file>>temp;file>>temp;file>>temp;

	file>>temp;
	getline(file,temp);

}

void extract_FeatureVector(char*fileName,int& n,int dim,double**& featureArray,double*& outputArray,bool isModelFile,SVM_stat& stat)
{
	fstream file;
	file.open(fileName);
	if(file.is_open()==false)
	{
		cout<<"Error in opening the file!"<<endl;
		return;
	}

	if(isModelFile==true)
	{
		extractModelStatistics(file,stat);
		n=stat.total_sv;
	}

	initArray(featureArray,outputArray,n,dim);

	string lineString;
	char*token;
	int d;
	double d_Value;
	for(int i=0;i<n;i++)
	{
		getline(file,lineString);

		//get the first token from the lineString
		token=strtok((char*)lineString.c_str()," :");
		outputArray[i]=atof(token);

		//Use tokenization to obtain the feature Array
		while(true)
		{
			token=strtok(NULL," :");

			//cout<<token<<endl;
			if(token!=NULL)
			{
				d=atoi(token);
				if(d==0)
					break;

				token=strtok(NULL," :");
				if(token!=NULL)
					d_Value=atof(token);
				else
				{
					cout<<"Read File ERROR!"<<endl;
					exit(1);
				}
				featureArray[i][d-1]=d_Value;
			}
			else
				break;
		}

	}

	file.close();

}

void outputResultFile(char*classResultName,SVM_stat& stat)
{
	fstream resultFile;
	resultFile.open(classResultName,ios::in|ios::out|ios::trunc);

	if(resultFile.is_open()==false)
	{
		cout<<"Cannot Open Result File!"<<endl;
		return;
	}

	for(int r=0;r<(int)stat.class_resultVector.size();r++)
		resultFile<<stat.class_resultVector[r]<<endl;

	resultFile.close();
}

void outputCheck(double** featureArray,double* outputArray,int n,int dim)
{
	for(int i=0;i<n;i++)
	{
		for(int d=0;d<dim;d++)
			cout<<featureArray[i][d]<<" ";
		cout<<endl;
	}

	cout<<endl;
	
	for(int i=0;i<n;i++)
		cout<<outputArray[i]<<" ";

	cout<<endl;
}

void saveMatrix_toFile(char*fileName,double**matrix,int numOfElements,int dim)
{
	fstream file;
	file.open(fileName,ios::in|ios::out|ios::trunc);
	if(file.is_open()==false)
	{
		cout<<"Cannot Open Files!"<<endl;
		return;
	}

	for(int i=0;i<numOfElements;i++)
	{
		for(int d=0;d<dim;d++)
			file<<matrix[i][d]<<" ";
		file<<endl;
	}

	file.close();
}

#ifdef EXACT_VALUE_STATS
void output_exact_ResultFile(char*classResultName,vector<double>& exactValueVector)
{
	fstream resultFile;
	resultFile.open(classResultName,ios::in|ios::out|ios::trunc);

	if(resultFile.is_open()==false)
	{
		cout<<"Cannot Open Result File!"<<endl;
		return;
	}

	for(int r=0;r<(int)exactValueVector.size();r++)
		resultFile<<exactValueVector[r]<<endl;

	resultFile.close();
}
#endif