#include "kd_tree.h"

double kdNode::LB(double*q,int dim,SVM_stat& stat)
{
	double lb,ub;
	double L_pos;
	double U_neg;

	lb=ell_MBR(q,boundary_neg,dim);
	ub=u_MBR(q,boundary_pos,dim);

	L_pos=sum_alpha_pos*exp(-stat.gammaValue*ub*ub); //LB_{G_{+}}
	U_neg=sum_alpha_neg*exp(-stat.gammaValue*lb*lb); //UB_{G_{-}}

	return (L_pos-U_neg);
}

double kdNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb,ub;
	double U_pos;
	double L_neg;

	lb=ell_MBR(q,boundary_pos,dim);
	ub=u_MBR(q,boundary_neg,dim);

	U_pos=sum_alpha_pos*exp(-stat.gammaValue*lb*lb);
	L_neg=sum_alpha_neg*exp(-stat.gammaValue*ub*ub);

	return (U_pos-L_neg);
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

//kdLinearAugNode
double kdLinearAugNode::LB(double*q,int dim,SVM_stat& stat)
{
	//LB_pos 
	double LB_pos;
	double t_star;
	double ip=0;

	//UB_neg
	double UB_neg;
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//obtain gamma_sum_pos and gamma_sum_neg
	//****************************************************************************//
	for(int d=0;d<dim;d++)
		ip=ip+q[d]*a_G_pos[d];
	gamma_sum_pos=stat.gammaValue*(sum_alpha_pos*stat.qSquareNorm-2*ip+S_G_pos);

	ip=0;
	for(int d=0;d<dim;d++)
		ip=ip+q[d]*a_G_neg[d];
	gamma_sum_neg=stat.gammaValue*(sum_alpha_neg*stat.qSquareNorm-2*ip+S_G_neg);
	//****************************************************************************//

	//Compute t^*
	if(sum_alpha_pos<epsilon)
		LB_pos=0;
	else
	{
		t_star=(gamma_sum_pos/sum_alpha_pos);
		//Compute LB_pos
		LB_pos=sum_alpha_pos*exp(-t_star);
	}

	if(sum_alpha_neg<epsilon)
		UB_neg=0;
	else
	{
		//Compute UB_neg
		//l_tri and u_tri
		/*ub=u_tri(q,O_r_neg,dim,radius_neg,temp_obt_dist_neg);
		if(temp_obt_dist_neg<radius_neg)
			lb=0;
		else
			lb=temp_obt_dist_neg-radius_neg;*/

		//l_MBR and u_MBR
		lb=ell_MBR(q,boundary_neg,dim);
		ub=u_MBR(q,boundary_neg,dim);

		l2=lb*lb;
		u2=ub*ub;

		//l2=u2, we cannot use mx+c to act as upper bound computationally
		if(u2-l2<epsilon)
			UB_neg=sum_alpha_neg*exp(-stat.gammaValue*l2);
		else
		{
			exp_l2=exp(-stat.gammaValue*l2);
			exp_u2=exp(-stat.gammaValue*u2);
			//compute m and c
			m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
			c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

			UB_neg=m*gamma_sum_neg+c*sum_alpha_neg;
		}
	}
	
	return (LB_pos-UB_neg);
}

double kdLinearAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	//LB_neg
	double LB_neg;
	double t_star;
	double ip=0;

	//UB_pos
	double UB_pos;
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//Compute t^*
	if(sum_alpha_neg<epsilon)
		LB_neg=0;
	else
	{
		t_star=(gamma_sum_neg/sum_alpha_neg);
		//Compute LB_neg
		LB_neg=sum_alpha_neg*exp(-t_star);
	}
	
	//Compute UB_pos
	if(sum_alpha_pos<epsilon)
		UB_pos=0;
	else
	{
		//l_tri and u_tri
		/*ub=u_tri(q,O_r_pos,dim,radius_pos,temp_obt_dist_pos);
		if(temp_obt_dist_pos<radius_pos)
			lb=0;
		else
			lb=temp_obt_dist_pos-radius_pos;*/

		//l_MBR and u_MBR
		lb=ell_MBR(q,boundary_pos,dim);
		ub=u_MBR(q,boundary_pos,dim);

		l2=lb*lb;
		u2=ub*ub;

		//l2=u2, we cannot use mx+c to act as upper bound computationally
		if(u2-l2<epsilon)
			UB_pos=sum_alpha_pos*exp(-stat.gammaValue*l2);
		else
		{
			exp_l2=exp(-stat.gammaValue*l2);
			exp_u2=exp(-stat.gammaValue*u2);
			//compute m and c
			m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
			c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

			UB_pos=m*gamma_sum_pos+c*sum_alpha_pos;
		}
	}
	
	return (UB_pos-LB_neg);
}

kdLinearAugNode*kdLinearAugNode::createNode()
{
	return new kdLinearAugNode();
}

void kdLinearAugNode::update_Aug(Node*node,Tree*t)
{
	this->update_linearAugInfo((kdLinearAugNode*)node,t);
}

void kdLinearAugNode::update_linearAugInfo(kdLinearAugNode*node,Tree*t)
{
	int id;
	double norm_square;
	//obtain vec a_G
	node->a_G_pos=new double[t->dim];
	node->a_G_neg=new double[t->dim];

	for(int d=0;d<t->dim;d++)
	{
		node->a_G_pos[d]=0;
		node->a_G_neg[d]=0;
	}

	//obtain a_G S_G (pos,neg)
	node->S_G_pos=0;
	node->S_G_neg=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		norm_square=0;
		for(int d=0;d<t->dim;d++)
			norm_square=norm_square+t->dataMatrix[id][d]*t->dataMatrix[id][d];

		if(t->outputArray[id]>0)
		{
			node->S_G_pos+=t->outputArray[id]*norm_square;
			for(int d=0;d<t->dim;d++)
				node->a_G_pos[d]+=t->outputArray[id]*t->dataMatrix[id][d];
		}
		else
		{
			node->S_G_neg+=fabs(t->outputArray[id])*norm_square;
			for(int d=0;d<t->dim;d++)
				node->a_G_neg[d]+=fabs(t->outputArray[id])*t->dataMatrix[id][d];
		}
	}

	if((int)node->idList.size()<=t->leafCapacity) //this is the leaf node
		return;

	update_linearAugInfo((kdLinearAugNode*)node->childVector[0],t);
	update_linearAugInfo((kdLinearAugNode*)node->childVector[1],t);
}

//kd-tree
kdTree::kdTree(int dim,double**dataMatrix,double*outputArray,int leafCapacity,SVM_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->leafCapacity=leafCapacity;
	this->outputArray=outputArray;
	this->stat=stat;
}

void kdTree::getNode_Boundary(kdNode*node)
{
	//double sum_alpha;
	int id;
	node->boundary_pos=new double*[dim];
	node->boundary_neg=new double*[dim];
	for(int d=0;d<dim;d++)
	{
		node->boundary_pos[d]=new double[2];
		node->boundary_neg[d]=new double[2];
	}

	for(int d=0;d<dim;d++)
	{
		node->boundary_pos[d][0]=inf;
		node->boundary_pos[d][1]=-inf;

		node->boundary_neg[d][0]=inf;
		node->boundary_neg[d][1]=-inf;
	}

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(outputArray[id]>=0)
		{
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[id][d]<node->boundary_pos[d][0])
					node->boundary_pos[d][0]=dataMatrix[id][d];
				if(dataMatrix[id][d]>node->boundary_pos[d][1])
					node->boundary_pos[d][1]=dataMatrix[id][d];
			}
		}
		else
		{
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[id][d]<node->boundary_neg[d][0])
					node->boundary_neg[d][0]=dataMatrix[id][d];
				if(dataMatrix[id][d]>node->boundary_neg[d][1])
					node->boundary_neg[d][1]=dataMatrix[id][d];
			}
		}
	}
}

double kdTree::obtain_SplitValue(kdNode*node,int split_Dim)
{
	vector<double> tempVector;
	int id;
	int middle_left,middle_right,middle;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		tempVector.push_back(dataMatrix[id][split_Dim]);
	}

	sort(tempVector.begin(),tempVector.end());

	if((int)tempVector.size()%2==0)//even number
	{
		middle_right=(int)tempVector.size()/2;
		middle_left=middle_right-1;

		return ((tempVector[middle_left]+tempVector[middle_right])/2.0);
	}
	else
	{
		middle=((int)tempVector.size()-1)/2;
		return tempVector[middle];
	}

	tempVector.clear();
}

void kdTree::KD_Tree_Recur(kdNode*node,int split_Dim)
{
	int id;
	int counter;
	//base case
	if((int)node->idList.size()<=leafCapacity)
	{
		node->initModel_inMemory(dataMatrix,outputArray,dim,stat);
		node->createModel_inMemory(dataMatrix,outputArray,dim,stat);
		return;
	}
	
	//added for testing
	//split_Dim=select_split_Dim(node);
	double splitValue=obtain_SplitValue(node,split_Dim); //code here

	//create two children
	kdNode*leftNode;
	kdNode*rightNode;

	leftNode=node->createNode();
	rightNode=node->createNode();

	counter=0;
	int halfSize=((int)node->idList.size())/2;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(dataMatrix[id][split_Dim]<=splitValue && counter<=halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	getNode_Boundary(leftNode);
	getNode_Boundary(rightNode);

	KD_Tree_Recur(leftNode,(split_Dim+1)%dim);
	KD_Tree_Recur(rightNode,(split_Dim+1)%dim);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree::build_kdTree(SVM_stat& stat)
{
	for(int i=0;i<stat.total_sv;i++)
		rootNode->idList.push_back(i);

	getNode_Boundary((kdNode*)rootNode);
	KD_Tree_Recur((kdNode*)rootNode,0);
	initTree_alpha((kdNode*)rootNode);
}

void kdTree::initTree_alpha(kdNode*node)
{
	double sum_alpha_pos=0;
	double sum_alpha_neg=0;
	int id;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(outputArray[id]>0)
			sum_alpha_pos+=outputArray[id];
		else
			sum_alpha_neg+=fabs(outputArray[id]);
	}

	node->sum_alpha_pos=sum_alpha_pos;
	node->sum_alpha_neg=sum_alpha_neg;

	if((int)node->idList.size()>leafCapacity)
	{
		initTree_alpha((kdNode*)node->childVector[0]);
		initTree_alpha((kdNode*)node->childVector[1]);
	}
}

void kdTree::updateAugment(kdNode*node)
{
	node->update_Aug(node,this);
}

//kd-tree_adv
void kdTree_adv::KD_Tree_adv_Recur(kdNode*node,int split_Dim)
{
	int id;
	int counter;
	//base case
	if((int)node->idList.size()<=leafCapacity)
	{
		getNode_Boundary(node);
		node->initModel_inMemory(dataMatrix,outputArray,dim,stat);
		node->createModel_inMemory(dataMatrix,outputArray,dim,stat);
		return;
	}
	
	//added for testing
	double splitValue=obtain_SplitValue(node,split_Dim); //code here

	//create two children
	kdNode*leftNode;
	kdNode*rightNode;

	leftNode=node->createNode();
	rightNode=node->createNode();

	counter=0;
	int halfSize=((int)node->idList.size())/2;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(dataMatrix[id][split_Dim]<=splitValue && counter<=halfSize)
		{
			leftNode->idList.push_back(id);
			counter++;
		}
		else
			rightNode->idList.push_back(id);
	}

	KD_Tree_adv_Recur(leftNode,(split_Dim+1)%dim);
	KD_Tree_adv_Recur(rightNode,(split_Dim+1)%dim);

	//code here
	merge_Boundary(node,leftNode->boundary_pos,rightNode->boundary_pos,true);
	merge_Boundary(node,leftNode->boundary_neg,rightNode->boundary_neg,false);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree_adv::build_kdTree_adv(SVM_stat& stat)
{
	for(int i=0;i<stat.total_sv;i++)
		rootNode->idList.push_back(i);

	KD_Tree_adv_Recur((kdNode*)rootNode,0);
	initTree_alpha((kdNode*)rootNode);
}

void kdTree_adv::merge_Boundary(kdNode*node,double**boundary_1,double**boundary_2,bool isPos)
{
	if(isPos==true)
	{
		node->boundary_pos=new double*[dim];
		for(int d=0;d<dim;d++)
			node->boundary_pos[d]=new double[2];

		//merge two boundaries
		for(int d=0;d<dim;d++)
		{
			if(boundary_1[d][0]<boundary_2[d][0])
				node->boundary_pos[d][0]=boundary_1[d][0];
			else
				node->boundary_pos[d][0]=boundary_2[d][0];

			if(boundary_1[d][1]<boundary_2[d][1])
				node->boundary_pos[d][1]=boundary_2[d][1];
			else
				node->boundary_pos[d][1]=boundary_1[d][1];
		}
	}
	else
	{
		node->boundary_neg=new double*[dim];
		for(int d=0;d<dim;d++)
			node->boundary_neg[d]=new double[2];

		//merge two boundaries
		for(int d=0;d<dim;d++)
		{
			if(boundary_1[d][0]<boundary_2[d][0])
				node->boundary_neg[d][0]=boundary_1[d][0];
			else
				node->boundary_neg[d][0]=boundary_2[d][0];

			if(boundary_1[d][1]<boundary_2[d][1])
				node->boundary_neg[d][1]=boundary_2[d][1];
			else
				node->boundary_neg[d][1]=boundary_1[d][1];
		}
	}
}

//kdLinearAugNode_adv
kdLinearAugNode_adv*kdLinearAugNode_adv::createNode()
{
	return new kdLinearAugNode_adv();
}

void kdLinearAugNode_adv::update_linearAugInfo_adv(kdLinearAugNode*node,Tree*t)
{
	int id;
	double norm_square_pos;
	double norm_square_neg;
	int idList_size=(int)node->idList.size();

	//obtain vec a_G
	node->a_G_pos=new double[t->dim];
	node->a_G_neg=new double[t->dim];

	if(idList_size<=t->leafCapacity)
	{
		node->S_G_pos=0;
		node->S_G_neg=0;

		for(int d=0;d<t->dim;d++)
		{
			node->a_G_pos[d]=0;
			node->a_G_neg[d]=0;
		}

		for(int i=0;i<idList_size;i++)
		{
			id=node->idList[i];

			if(t->outputArray[id]>0)
			{
				for(int d=0;d<t->dim;d++)
					node->a_G_pos[d]+=t->outputArray[id]*t->dataMatrix[id][d];

				norm_square_pos=0;
				for(int d=0;d<t->dim;d++)
					norm_square_pos=norm_square_pos+t->dataMatrix[id][d]*t->dataMatrix[id][d];

				node->S_G_pos+=t->outputArray[id]*norm_square_pos;
			}
			else
			{
				for(int d=0;d<t->dim;d++)
					node->a_G_neg[d]+=fabs(t->outputArray[id])*t->dataMatrix[id][d];

				norm_square_neg=0;
				for(int d=0;d<t->dim;d++)
					norm_square_neg=norm_square_neg+t->dataMatrix[id][d]*t->dataMatrix[id][d];

				node->S_G_neg+=fabs(t->outputArray[id])*norm_square_neg;
			}
		}
		return;
	}
	
	update_linearAugInfo_adv((kdLinearAugNode_adv*)node->childVector[0],t);
	update_linearAugInfo_adv((kdLinearAugNode_adv*)node->childVector[1],t);

	//update a_G and S_G in the internal nodes
	for(int d=0;d<t->dim;d++)
	{
		node->a_G_pos[d]=((kdLinearAugNode_adv*)node->childVector[0])->a_G_pos[d]+((kdLinearAugNode_adv*)node->childVector[1])->a_G_pos[d];
		node->a_G_neg[d]=((kdLinearAugNode_adv*)node->childVector[0])->a_G_neg[d]+((kdLinearAugNode_adv*)node->childVector[1])->a_G_neg[d];
	}
	
	node->S_G_pos=((kdLinearAugNode_adv*)node->childVector[0])->S_G_pos+((kdLinearAugNode_adv*)node->childVector[1])->S_G_pos;
	node->S_G_neg=((kdLinearAugNode_adv*)node->childVector[0])->S_G_neg+((kdLinearAugNode_adv*)node->childVector[1])->S_G_neg;
}

void kdLinearAugNode_adv::update_Aug(Node*node,Tree*t)
{
	this->update_linearAugInfo_adv((kdLinearAugNode*)node,t);
}