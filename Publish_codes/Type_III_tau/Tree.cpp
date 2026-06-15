#include "Tree.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void Node::initModel_inMemory(double**dataMatrix,double*outputArray,int dim,SVM_stat& stat)
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
	model->l=(int)idList.size();

	int elements=0;

	for(int i=0;i<(int)idList.size();i++)
	{
		for(int d=0;d<dim;d++)
		{
			if(fabs(dataMatrix[idList[i]][d])>epsilon)
				elements++;
		}
	}
		
	elements+=(int)idList.size();

	x_space = NULL;
	x_space=Malloc(svm_node,elements);

	model->sv_coef = Malloc(double *,1);
	model->sv_coef[0] = Malloc(double,(int)idList.size());
	model->SV = Malloc(svm_node*,(int)idList.size());
}

void Node::createModel_inMemory(double**dataMatrix,double*outputArray,int dim,SVM_stat& stat)
{
	int pos=0;
	for(int i=0;i<(int)idList.size();i++)
	{
		model->sv_coef[0][i]=outputArray[idList[i]];
		model->SV[i]=&x_space[pos];
		for(int d=0;d<dim;d++)
		{
			if(fabs(dataMatrix[idList[i]][d])>epsilon)
			{
				x_space[pos].index=d+1;
				x_space[pos].value=dataMatrix[idList[i]][d];
				pos++;
			}
		}
		x_space[pos].index=-1;
		pos++;
	}
}

void Node::freeModel_inMemory()
{
	svm_free_and_destroy_model(&this->model);
}