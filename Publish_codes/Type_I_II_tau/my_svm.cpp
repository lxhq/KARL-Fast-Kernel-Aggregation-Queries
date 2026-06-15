#include "my_svm.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void convert_2D_qMatrix_to_SVMFormat(double**queryMatrix,svm_node**& svm_qMatrix,int dim,int qNum)
{
	int position;
	svm_qMatrix=new svm_node*[qNum];
	for(int q=0;q<qNum;q++)
		svm_qMatrix[q]=new svm_node[dim];

	for(int q=0;q<qNum;q++)
	{
		position=0;
		for(int d=0;d<dim;d++)
		{
			if(fabs(queryMatrix[q][d])>epsilon)
			{
				svm_qMatrix[q][position].index=d+1;
				svm_qMatrix[q][position].value=queryMatrix[q][d];
				position++;
			}
		}
		svm_qMatrix[q][position].index=-1;
	}
}

void loadSVM_model(char*modelName,svm_model*& model)
{
	if((model=svm_load_model(modelName))==0)
	{
		fprintf(stderr,"can't open model file %s\n","");
		exit(1);
	}
}

void reorder_model(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model)
{
	fstream model_file;

	model_file.open("reorder_model.txt",ios::in | ios::out | ios::trunc);
	if(model_file.is_open()==false)
	{
		cout<<"Cannot store new (re-ordered) model file!"<<endl;
		exit(1);
	}
	//store to model_file
	//**********************************************//
	model_file<<"svm_type one_class"<<endl;
	model_file<<"kernel_type rbf"<<endl;
	model_file<<"gamma "<<stat.gammaValue<<endl;
	model_file<<"nr_class 2"<<endl;
	model_file<<"total_sv "<<stat.total_sv<<endl;
	model_file<<"rho "<<stat.rho<<endl;
	model_file<<"SV"<<endl;

	for(int i=0;i<stat.total_sv;i++)
	{
		model_file<<alphaArray[i]<<" ";
		for(int d=0;d<dim;d++)
		{
			if(P[i][d]>epsilon)
				model_file<<(d+1)<<":"<<P[i][d]<<" ";
		}
		model_file<<endl;
	}

	model_file.close();
	//**********************************************//

	loadSVM_model((char*)"reorder_model.txt",model);
}

void init_model_inMemory(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space)
{
	model=Malloc(svm_model,1);
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->sv_indices = NULL;
	model->label = NULL;
	model->nSV = NULL;

	//update the parameters in svm model
	model->param.svm_type=2;
	model->param.kernel_type=2;
	model->param.gamma=stat.gammaValue;
	model->nr_class=2;
	model->rho=new double[1];
	model->rho[0]=stat.rho;
	model->l=stat.total_sv;

	int elements=0;
	for(int i=0;i<stat.total_sv;i++)
		for(int d=0;d<dim;d++)
			if(fabs(P[i][d])>epsilon)
				elements++;
	elements+=stat.total_sv;

	x_space = NULL;
	x_space=Malloc(svm_node,elements);

	model->sv_coef = Malloc(double *,1);
	model->sv_coef[0] = Malloc(double,stat.total_sv);
	model->SV = Malloc(svm_node*,stat.total_sv);
}

void createModel_inMemory(double**P,double*alphaArray,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space)
{
	int pos=0;
	for(int i=0;i<stat.total_sv;i++)
	{
		model->sv_coef[0][i]=alphaArray[i];
		model->SV[i]=&x_space[pos];
		for(int d=0;d<dim;d++)
		{
			if(fabs(P[i][d])>epsilon)
			{
				x_space[pos].index=d+1;
				x_space[pos].value=P[i][d];
				pos++;
			}
		}
		x_space[pos].index=-1;
		pos++;
	}
}

void prediction(const svm_model *model, const svm_node *x,SVM_stat& stat)
{
	double predictValue=svm_predict(model,x);
	stat.class_resultVector.push_back((int)predictValue);
}