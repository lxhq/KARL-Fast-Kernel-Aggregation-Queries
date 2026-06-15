#include "binaryTree.h"
#include "GBF_KAQ.h"

struct order_Entry
{
	int id;
	double dist;
	double norm;
	double ip;
};

bool compare_order_Entry(const order_Entry& o_Entry1,const order_Entry& o_Entry2)
{
	return o_Entry1.dist < o_Entry2.dist; 
}

bool compare_order_Entry_decr(const order_Entry& o_Entry1,const order_Entry& o_Entry2)
{
	return o_Entry1.dist > o_Entry2.dist;
}

bool compare_order_Entry_norm(const order_Entry& o_Entry1,const order_Entry& o_Entry2)
{
	return o_Entry1.norm < o_Entry2.norm;
}

bool compare_order_Entry_ip(const order_Entry& o_Entry1,const order_Entry& o_Entry2)
{
	return o_Entry1.ip > o_Entry2.ip;
}

//binaryNode
double binaryNode::LB(double*q,int dim,KDE_stat& stat)
{
	double t_star;
	double ip=0;
	for(int d=0;d<dim;d++)
		ip=ip+q[d]*a_G[d];

	//Compute t^*
	gamma_sum=sumE*stat.qSquareNorm-2*ip+S_G;
	t_star=gamma_sum/sumE;

	return sumE*exp(-t_star);
}

double binaryNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//l_tri and u_tri
	ub=u_tri(q,center,dim,radius,temp_obt_dist);
	if(temp_obt_dist<radius)
		lb=0;
	else
		lb=temp_obt_dist-radius;

	l2=lb*lb;
	u2=ub*ub;

	//l2=u2, we cannot use mx+c to act as upper bound computationally
	if(u2-l2<epsilon)
		return sumE*exp(-l2);

	exp_l2=exp(-l2);
	exp_u2=exp(-u2);
	//compute m and c
	m=(exp_u2-exp_l2)/(u2-l2);
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*sumE);
}

binaryNode*binaryNode::createNode()
{
	return new binaryNode();
}

void binaryNode::update_Aug(Node*node,Tree*t)
{
	//no code
}

void binaryNode::update_Aug(double**dataMatrix,int dim,KDE_stat& stat)
{
	int id;
	double temp_norm_sq;
	double temp_radius;
	double max_radius=-inf;
	//initialization 
	sumE=0;
	a_G=new double[dim];
	center=new double[dim];
	for(int d=0;d<dim;d++)
	{
		a_G[d]=0;
		center[d]=0;
	}
	S_G=0;

	//(1)update sumE //(2)update a_G //(3)update S_G //(4)update center //(5)update radius
	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		sumE+=1; //(1)
		temp_norm_sq=0;
		for(int d=0;d<dim;d++)
		{
			a_G[d]+=dataMatrix[id][d]; //(2)
			temp_norm_sq+=dataMatrix[id][d]*dataMatrix[id][d];
		}
		S_G+=temp_norm_sq; //(3)
	}

	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sumE; //(4)

	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		temp_radius=0;
		for(int d=0;d<dim;d++)
			temp_radius+=(dataMatrix[id][d]-center[d])*(dataMatrix[id][d]-center[d]);

		if(temp_radius>max_radius)
			max_radius=temp_radius;
	}

	radius=sqrt(max_radius); //(5)
}

void binaryNode::Aug_Incr(double**dataMatrix,int dim,KDE_stat& stat,int id)
{
	double square_norm=0;
	
	//(1)update sumE //(2)update a_G //(3)update S_G //(4)update center //(5)update radius
	sumE+=1; //(1)
	for(int d=0;d<dim;d++)
	{
		a_G[d]+=dataMatrix[id][d]; //(2)
		square_norm+=dataMatrix[id][d]*dataMatrix[id][d]; 
	}
	S_G+=square_norm; //(3)
	
	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sumE; //(4)
}

void binaryNode::Aug_Decr(double**dataMatrix,int dim,KDE_stat& stat,int id)
{
	double square_norm=0;
	//(1)update sumE //(2)update a_G //(3)update S_G //(4)update center
	sumE-=1; //(1)
	for(int d=0;d<dim;d++)
	{
		a_G[d]-=dataMatrix[id][d]; //(2)
		square_norm+=dataMatrix[id][d]*dataMatrix[id][d]; 
	}

	S_G-=square_norm; //(3)
	
	if((int)idList.size()>0)
	{
		for(int d=0;d<dim;d++)
			center[d]=a_G[d]/sumE; //(4)
	}
	else
	{
		for(int d=0;d<dim;d++)
			center[d]=0; //(4)
	}
}

void binaryNode::update_Radius(double**dataMatrix,int dim)
{
	int id;
	double temp_radius;
	double max_radius=-inf;

	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		temp_radius=0;
		for(int d=0;d<dim;d++)
			temp_radius+=(dataMatrix[id][d]-center[d])*(dataMatrix[id][d]-center[d]);

		if(temp_radius>max_radius)
			max_radius=temp_radius;
	}

	radius=sqrt(max_radius);
}

void half_division(vector<order_Entry>& o_EntryVector,binaryNode*node1,binaryNode*node2,double**dataMatrix,int dim,KDE_stat& stat)
{
	int halfSize=(int)ceil((double)o_EntryVector.size()/2.0);
	for(int i=0;i<(int)o_EntryVector.size();i++)
	{
		if(i<halfSize)
		{
			node1->idList.push_back(o_EntryVector[i].id);
			node1->Aug_Incr(dataMatrix,dim,stat,o_EntryVector[i].id);
		}
		else
		{
			node2->idList.push_back(o_EntryVector[i].id);
			node2->Aug_Incr(dataMatrix,dim,stat,o_EntryVector[i].id);
		}
	}

	node1->update_Radius(dataMatrix,dim);
	node2->update_Radius(dataMatrix,dim);

	o_EntryVector.clear();
}

void binaryNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,int dim,KDE_stat& stat)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		o_Entry.id=this->idList[i];
		o_Entry.dist=euclid_dist(dataMatrix[o_Entry.id],this->center,dim);
		o_EntryVector.push_back(o_Entry);
	}

	sort(o_EntryVector.begin(),o_EntryVector.end(),compare_order_Entry);

	half_division(o_EntryVector,node1,node2,dataMatrix,dim,stat);
}

void binaryNode::initNode(int dim)
{
	//(1)init sum_alpha //(2)init a_G //(3)init S_G //(4)init center //(5) radius
	sumE=0; //(1)
	S_G=0; //(3)
	radius=0; //(5)
	a_G=new double[dim];
	center=new double[dim];
	for(int d=0;d<dim;d++)
	{
		a_G[d]=0; //(2)
		center[d]=0; //(4)
	}
}

void binaryNode::assign(binaryNode*bNode_delta,int dim)
{
	//(1)assign sumE //(2)assign a_G //(3)assign S_G //(4)assign center //(5)assign radius
	sumE=bNode_delta->sumE; //(1)
	S_G=bNode_delta->S_G; //(3)
	//radius=bNode_delta->radius; //(5)
	for(int d=0;d<dim;d++)
	{
		a_G[d]=bNode_delta->a_G[d]; //(2)
		center[d]=bNode_delta->center[d]; //(4)
	}
}

//ballNode
ballNode*ballNode::createNode()
{
	return new ballNode();
}

double ballNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//l_tri and u_tri
	/*ub=u_tri(q,center,dim,radius,temp_obt_dist);
	if(temp_obt_dist<radius)
		lb=0;
	else
		lb=temp_obt_dist-radius;*/

	//l_MBR and u_MBR
	lb=ell_MBR(q,boundary,dim);
	ub=u_MBR(q,boundary,dim);

	l2=lb*lb;
	u2=ub*ub;

	//l2=u2, we cannot use mx+c to act as upper bound computationally
	if(u2-l2<epsilon)
		return sumE*exp(-l2);

	exp_l2=exp(-l2);
	exp_u2=exp(-u2);
	//compute m and c
	m=(exp_u2-exp_l2)/(u2-l2);
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*sumE);
}

void ballNode::updateBoundary(double**dataMatrix,KDE_stat& stat,int dim)
{
	int id;

	boundary=new double*[dim];
	for(int d=0;d<dim;d++)
		boundary[d]=new double[2];

	//initialization of boundary
	for(int d=0;d<dim;d++)
	{
		boundary[d][0]=inf;
		boundary[d][1]=-inf;
	}

	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		for(int d=0;d<dim;d++)
		{
			if(dataMatrix[id][d]<boundary[d][0])
				boundary[d][0]=dataMatrix[id][d];
			if(dataMatrix[id][d]>boundary[d][1])
				boundary[d][1]=dataMatrix[id][d];
		}
	}
}

void ballNode::update_Aug(double**dataMatrix,int dim,KDE_stat& stat)
{
	binaryNode::update_Aug(dataMatrix,dim,stat);
	updateBoundary(dataMatrix,stat,dim);
}

void find_furthest(double*pivot,vector<int>& idList,double**dataMatrix,int dim,int& id1,int& id2)
{
	double max_dist=-inf;
	double cur_dist;
	int cur_id;

	for(int i=0;i<(int)idList.size();i++)
	{
		cur_id=idList[i];
		cur_dist=euclid_dist(pivot,dataMatrix[cur_id],dim);

		if(cur_dist>max_dist)
		{
			max_dist=cur_dist;
			id1=cur_id;
		}
	}

	max_dist=-inf;
	for(int i=0;i<(int)idList.size();i++)
	{
		cur_id=idList[i];
		if(cur_id==id1)
			continue;

		cur_dist=euclid_dist(dataMatrix[id1],dataMatrix[cur_id],dim);

		if(cur_dist>max_dist)
		{
			max_dist=cur_dist;
			id2=cur_id;
		}
	}
}

void balance_saparation_basic(double**dataMatrix,int dim,int id1,int id2,vector<order_Entry>& o_EntryVector,binaryNode*node)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector1;
	vector<order_Entry> o_EntryVector2;
	double euclid_1;
	double euclid_2;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		o_Entry.id=node->idList[i];

		if((int)o_EntryVector1.size()>(int)(node->idList.size()/2.0))
		{
			o_Entry.dist=euclid_dist(dataMatrix[id2],dataMatrix[node->idList[i]],dim);
			o_EntryVector2.push_back(o_Entry);
			continue;
		}

		if((int)o_EntryVector2.size()>(int)(node->idList.size()/2.0))
		{
			o_Entry.dist=euclid_dist(dataMatrix[id1],dataMatrix[node->idList[i]],dim);
			o_EntryVector1.push_back(o_Entry);
			continue;
		}

		euclid_1=euclid_dist(dataMatrix[id1],dataMatrix[node->idList[i]],dim);
		euclid_2=euclid_dist(dataMatrix[id2],dataMatrix[node->idList[i]],dim);
		if(euclid_1<euclid_2)
		{
			o_Entry.dist=euclid_1;
			o_EntryVector1.push_back(o_Entry);
		}
		else
		{
			o_Entry.dist=euclid_2;
			o_EntryVector2.push_back(o_Entry);
		}
	}

	sort(o_EntryVector1.begin(),o_EntryVector1.end(),compare_order_Entry);
	sort(o_EntryVector2.begin(),o_EntryVector2.end(),compare_order_Entry_decr);

	o_EntryVector.insert(o_EntryVector.end(),o_EntryVector1.begin(),o_EntryVector1.end());
	o_EntryVector.insert(o_EntryVector.end(),o_EntryVector2.begin(),o_EntryVector2.end());
	o_EntryVector1.clear();
	o_EntryVector2.clear();
}

void ballNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,int dim,KDE_stat& stat)
{
	static int isRoot=1;
	int id1,id2;
	vector<order_Entry> o_EntryVector;

	if(isRoot==1)
		isRoot=0;
	else
		updateBoundary(dataMatrix,stat,dim);

	find_furthest(this->center,this->idList,dataMatrix,dim,id1,id2);
	balance_saparation_basic(dataMatrix,dim,id1,id2,o_EntryVector,this);

	half_division(o_EntryVector,node1,node2,dataMatrix,dim,stat);
}

ballNode_SOTA*ballNode_SOTA::createNode()
{
	return new ballNode_SOTA();
}

double ballNode_SOTA::LB(double*q,int dim,KDE_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sumE*exp(-ub*ub);

	return L;
}

double ballNode_SOTA::UB(double*q,int dim,KDE_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sumE*exp(-lb*lb);

	return U;
}

void binaryTree::build_BinaryTreeRecur(binaryNode*node)
{
	if((int)node->idList.size()>leafCapacity)
	{
		binaryNode*node1=node->createNode();
		binaryNode*node2=node->createNode();

		node1->initNode(dim);
		node2->initNode(dim);
		
		node->GPA(node1,node2,dataMatrix,dim,stat);

		build_BinaryTreeRecur(node1);
		build_BinaryTreeRecur(node2);

		node->childVector.push_back(node1);
		node->childVector.push_back(node2);
	}
	else //leafNode
		node->update_Aug(dataMatrix,dim,stat);
}

void binaryTree::build_BinaryTree()
{
	binaryNode*& bNode=(binaryNode*&) rootNode;
	for(int i=0;i<stat.n;i++)
		bNode->idList.push_back(i);

	bNode->update_Aug(dataMatrix,dim,stat);

	build_BinaryTreeRecur((binaryNode*)rootNode);
}

binaryTree::binaryTree(int dim,double**dataMatrix,int leafCapacity,KDE_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->leafCapacity=leafCapacity;
	this->stat=stat;
}