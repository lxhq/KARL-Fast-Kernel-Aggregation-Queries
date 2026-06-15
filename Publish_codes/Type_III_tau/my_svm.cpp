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

void createModel_inMemory(double**P,double*outputArray,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space)
{
	int pos=0;
	for(int i=0;i<stat.total_sv;i++)
	{
		model->sv_coef[0][i]=outputArray[i];
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

void init_model_inMemory(double**P,int dim,SVM_stat& stat,svm_model*& model,svm_node*& x_space)
{
	model=Malloc(svm_model,1);
	model->rho = NULL;
	model->probA = NULL;
	model->probB = NULL;
	model->sv_indices = NULL;
	model->label = NULL;
	model->nSV = NULL;

	//update the parameters in svm model
	model->param.svm_type=0;
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

void prediction(const svm_model *model, const svm_node *x,SVM_stat& stat)
{
	double predictValue=svm_predict(model,x);
	stat.class_resultVector.push_back((int)predictValue);
}