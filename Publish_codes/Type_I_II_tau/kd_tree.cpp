#include "kd_tree.h"

double kdNode::LB(double*q,int dim,SVM_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sum_alpha*exp(-stat.gammaValue*ub*ub);

	return L;
}

double kdNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sum_alpha*exp(-stat.gammaValue*lb*lb);

	return U;
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

//kdAugNode
double kdAugNode::LB(double*q,int dim,SVM_stat& stat)
{
	double ub;
	double L;

	//u_tri
	ub=u_tri(q,center,dim,radius,temp_obt_dist);

	L=sum_alpha*exp(-stat.gammaValue*ub*ub);

	return L;
}

double kdAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double U;

	if(temp_obt_dist<radius)
		lb=0;
	else
		lb=temp_obt_dist-radius;

	U=sum_alpha*exp(-stat.gammaValue*lb*lb);

	return U;
}

void kdAugNode::update_Aug(Node*node,Tree*t)
{
	updateAugInfo((kdAugNode*)node,t);
}

kdAugNode*kdAugNode::createNode()
{
	return new kdAugNode();
}
void kdAugNode::updateAugInfo(kdAugNode*node,Tree*t)
{
	double max_Radius=-inf;
	double temp_Radius;
	int id;

	node->center=new double[t->dim];
	for(int d=0;d<t->dim;d++)
		node->center[d]=0;

	//find center
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		for(int d=0;d<t->dim;d++)
			node->center[d]+=t->alphaArray[id]*t->dataMatrix[id][d];
	}
	for(int d=0;d<t->dim;d++)
		node->center[d]=node->center[d]/(double)node->sum_alpha;

	//find radius
	for(int i=0;i<(int)node->idList.size();i++)
	{
		temp_Radius=0;
		id=node->idList[i];
		for(int d=0;d<t->dim;d++)
			temp_Radius=temp_Radius+(t->dataMatrix[id][d]-node->center[d])*(t->dataMatrix[id][d]-node->center[d]);

		if(temp_Radius>max_Radius)
			max_Radius=temp_Radius;
	}
	node->radius=sqrt(max_Radius);

	if((int)node->idList.size()<=t->leafCapacity) //this is the leaf node
		return;

	updateAugInfo((kdAugNode*)node->childVector[0],t);
	updateAugInfo((kdAugNode*)node->childVector[1],t);
}
//No LB for this node
double kdAugNode_MIN::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double U;

	lb=ell_MIN(q,boundary,d_SortObj,(int)idList.size()-1,dim);

	U=sum_alpha*exp(-stat.gammaValue*lb*lb);

	return U;
}

kdAugNode_MIN*kdAugNode_MIN::createNode()
{
	return new kdAugNode_MIN();
}

void kdAugNode_MIN::update_Aug(Node*node,Tree*t)
{
	updateAugMinInfo((kdAugNode_MIN*)node,t);
}

void kdAugNode_MIN::updateAugMinInfo(kdAugNode_MIN*node,Tree*t)
{
	int id;
	node->d_SortObj=new double*[t->dim];
	for(int d=0;d<t->dim;d++)
		node->d_SortObj[d]=new double[(int)node->idList.size()];

	for(int d=0;d<t->dim;d++)
	{
		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			node->d_SortObj[d][i]=t->dataMatrix[id][d];
		}

		sort(node->d_SortObj[d],node->d_SortObj[d]+(int)node->idList.size());
	}

	if((int)node->idList.size()<=t->leafCapacity) //this is the leaf node
		return;

	updateAugMinInfo((kdAugNode_MIN*)node->childVector[0],t);
	updateAugMinInfo((kdAugNode_MIN*)node->childVector[1],t);
}
//kdLinearAugNode
double kdLinearAugNode::LB(double*q,int dim,SVM_stat& stat)
{
	double t_star;
	double ip=0;
	for(int d=0;d<dim;d++)
		ip=ip+q[d]*a_G[d];

	//Compute t^*
	gamma_sum=stat.gammaValue*(sum_alpha*stat.qSquareNorm-2*ip+S_G);
	t_star=(gamma_sum/sum_alpha);

	return sum_alpha*exp(-t_star);
}

double kdLinearAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//use MBR
	lb=ell_MBR(q,boundary,dim);
	ub=u_MBR(q,boundary,dim);

	//l_tri and u_tri (don't use triangle inequality)
	/*ub=u_tri(q,center,dim,radius,temp_obt_dist);
	if(temp_obt_dist<radius)
		lb=0;
	else
		lb=temp_obt_dist-radius;*/

	l2=lb*lb;
	u2=ub*ub;

	//l2=u2, we cannot use mx+c to act as upper bound computationally
	if(u2-l2<epsilon)
		return sum_alpha*exp(-stat.gammaValue*l2);

	exp_l2=exp(-stat.gammaValue*l2);
	exp_u2=exp(-stat.gammaValue*u2);
	//compute m and c
	m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*sum_alpha);
}

kdLinearAugNode*kdLinearAugNode::createNode()
{
	return new kdLinearAugNode();
}

void kdLinearAugNode::update_Aug(Node*node,Tree*t)
{
	this->updateAugInfo((kdLinearAugNode*)node,t);
	this->update_linearAugInfo((kdLinearAugNode*)node,t);
}

void kdLinearAugNode::update_linearAugInfo(kdLinearAugNode*node,Tree*t)
{
	int id;
	double norm_square;
	//obtain vec a_G
	node->a_G=new double[t->dim];

	//We relate center and a_G in this way
	for(int d=0;d<t->dim;d++)
		node->a_G[d]=node->center[d]*node->sum_alpha;

	//obtain S_G
	node->S_G=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		norm_square=0;
		for(int d=0;d<t->dim;d++)
			norm_square=norm_square+t->dataMatrix[id][d]*t->dataMatrix[id][d];

		node->S_G+=t->alphaArray[id]*norm_square;
	}

	if((int)node->idList.size()<=t->leafCapacity) //this is the leaf node
		return;

	update_linearAugInfo((kdLinearAugNode*)node->childVector[0],t);
	update_linearAugInfo((kdLinearAugNode*)node->childVector[1],t);
}

//kd-tree
kdTree::kdTree(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->leafCapacity=leafCapacity;
	this->alphaArray=alphaArray;
	this->stat=stat;
}

void kdTree::getNode_Boundary(kdNode*node)
{
	//double sum_alpha;
	int id;
	node->boundary=new double*[dim];
	for(int d=0;d<dim;d++)
		node->boundary[d]=new double[2];

	for(int d=0;d<dim;d++)
	{
		node->boundary[d][0]=inf;
		node->boundary[d][1]=-inf;
	}

	for(int d=0;d<dim;d++)
	{
		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			if(dataMatrix[id][d]<node->boundary[d][0])
				node->boundary[d][0]=dataMatrix[id][d];

			if(dataMatrix[id][d]>node->boundary[d][1])
				node->boundary[d][1]=dataMatrix[id][d];
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
		node->initModel_inMemory(dataMatrix,alphaArray,dim,stat);
		node->createModel_inMemory(dataMatrix,alphaArray,dim,stat);
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

int kdTree::select_split_Dim(kdNode*node)
{
	double interval_difference=0;
	double largest_difference=0;
	int best_Dim=-1;

	for(int d=0;d<dim;d++)
	{
		interval_difference=(node->boundary[d][1]-node->boundary[d][0]);
		if(largest_difference<interval_difference)
		{
			largest_difference=interval_difference;
			best_Dim=d;
		}
	}

	return best_Dim;
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
	double sum_alpha=0;
	int id;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		sum_alpha+=alphaArray[id];
	}
	node->sum_alpha=sum_alpha;

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

void kdTree::lookUp_KDTree()
{
	lookUp_KDTree_Recur((kdNode*)rootNode);
}

void kdTree::lookUp_KDTree_Recur(kdNode*node)
{
	int idListSize=node->idList.size();
	for(int i=0;i<idListSize;i++)
		cout<<node->idList[i]<<" ";
	cout<<endl;
	cout<<"sum_alpha: "<<node->sum_alpha<<endl;

	if(idListSize>leafCapacity)
	{
		lookUp_KDTree_Recur((kdNode*)node->childVector[0]);
		lookUp_KDTree_Recur((kdNode*)node->childVector[1]);
	}
	else
	{
		int id;
		for(int i=0;i<idListSize;i++)
		{
			id=node->idList[i];
			cout<<"(";
			for(int d=0;d<dim-1;d++)
				cout<<dataMatrix[id][d]<<",";
			cout<<dataMatrix[id][dim-1]<<")"<<endl;
		}
	}
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
		node->initModel_inMemory(dataMatrix,alphaArray,dim,stat);
		node->createModel_inMemory(dataMatrix,alphaArray,dim,stat);

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

	//get back the boundary of leftNode and rightNode using O(d) time by merging boundaries from left and right node
	merge_Boundary(node,leftNode->boundary,rightNode->boundary);

	node->childVector.push_back(leftNode);
	node->childVector.push_back(rightNode);
}

void kdTree_adv::merge_Boundary(kdNode*node,double**boundary_1,double**boundary_2)
{
	node->boundary=new double*[dim];
	for(int d=0;d<dim;d++)
		node->boundary[d]=new double[2];

	//merge two boundaries
	for(int d=0;d<dim;d++)
	{
		if(boundary_1[d][0]<boundary_2[d][0])
			node->boundary[d][0]=boundary_1[d][0];
		else
			node->boundary[d][0]=boundary_2[d][0];

		if(boundary_1[d][1]<boundary_2[d][1])
			node->boundary[d][1]=boundary_2[d][1];
		else
			node->boundary[d][1]=boundary_1[d][1];
	}
}

void kdTree_adv::build_kdTree_adv(SVM_stat& stat)
{
	for(int i=0;i<stat.total_sv;i++)
		rootNode->idList.push_back(i);

	KD_Tree_adv_Recur((kdNode*)rootNode,0);
	initTree_alpha((kdNode*)rootNode);
}

//kdLinearAugNode_adv
kdLinearAugNode_adv*kdLinearAugNode_adv::createNode()
{
	return new kdLinearAugNode_adv();
}

void kdLinearAugNode_adv::update_linearAugInfo_adv(kdLinearAugNode*node,Tree*t)
{
	int id;
	double norm_square;
	int idList_size=(int)node->idList.size();

	//obtain vec a_G
	node->a_G=new double[t->dim];

	if(idList_size<=t->leafCapacity)
	{
		node->S_G=0;
		for(int d=0;d<t->dim;d++)
			node->a_G[d]=0;

		for(int i=0;i<idList_size;i++)
		{
			id=node->idList[i];
			for(int d=0;d<t->dim;d++)
				node->a_G[d]+=t->alphaArray[id]*t->dataMatrix[id][d];

			norm_square=0;
			for(int d=0;d<t->dim;d++)
				norm_square=norm_square+t->dataMatrix[id][d]*t->dataMatrix[id][d];

			node->S_G+=t->alphaArray[id]*norm_square;
		}
		return;
	}
	
	update_linearAugInfo_adv((kdLinearAugNode_adv*)node->childVector[0],t);
	update_linearAugInfo_adv((kdLinearAugNode_adv*)node->childVector[1],t);

	//update a_G and S_G in the internal nodes
	for(int d=0;d<t->dim;d++)
		node->a_G[d]=((kdLinearAugNode_adv*)node->childVector[0])->a_G[d]+((kdLinearAugNode_adv*)node->childVector[1])->a_G[d];

	node->S_G=((kdLinearAugNode_adv*)node->childVector[0])->S_G+((kdLinearAugNode_adv*)node->childVector[1])->S_G;
}

void kdLinearAugNode_adv::update_Aug(Node*node,Tree*t)
{
	this->update_linearAugInfo_adv((kdLinearAugNode*)node,t);
}