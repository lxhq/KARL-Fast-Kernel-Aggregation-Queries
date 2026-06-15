#include "kd_tree.h"

double kdNode::LB(double*q,int dim,KDE_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sumE*exp(-ub*ub);

	return L;
}

double kdNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sumE*exp(-lb*lb);

	return U;
}

kdNode*kdNode::createNode()
{
	return new kdNode();
}

//kdAugNode
double kdAugNode::LB(double*q,int dim,KDE_stat& stat)
{
	double ub;
	double L;

	//u_tri
	ub=u_tri(q,center,dim,radius,temp_obt_dist);

	L=sumE*exp(-ub*ub);

	return L;
}

double kdAugNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb;
	double U;

	if(temp_obt_dist<radius)
		lb=0;
	else
		lb=temp_obt_dist-radius;

	U=sumE*exp(-lb*lb);

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
			node->center[d]+=t->dataMatrix[id][d];
	}
	for(int d=0;d<t->dim;d++)
		node->center[d]=node->center[d]/node->sumE;

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

//kdLinearAugNode
double kdLinearAugNode::LB(double*q,int dim,KDE_stat& stat)
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

double kdLinearAugNode::UB(double*q,int dim,KDE_stat& stat)
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
		node->a_G[d]=node->center[d]*node->sumE;

	//obtain S_G
	node->S_G=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		norm_square=0;
		for(int d=0;d<t->dim;d++)
			norm_square=norm_square+t->dataMatrix[id][d]*t->dataMatrix[id][d];

		node->S_G+=norm_square;
	}

	if((int)node->idList.size()<=t->leafCapacity) //this is the leaf node
		return;

	update_linearAugInfo((kdLinearAugNode*)node->childVector[0],t);
	update_linearAugInfo((kdLinearAugNode*)node->childVector[1],t);
}

static double expFromLogValue(long double log_value)
{
	if(log_value==numeric_limits<long double>::infinity())
		return inf;
	if(log_value==(-numeric_limits<long double>::infinity()))
		return 0.0;
	if(!isfinite((double)log_value))
		return inf;
	if(log_value<-745.0L)
		return 0.0;
	if(log_value>700.0L)
		return inf;

	return (double)expl(log_value);
}

static double logSumExp2(long double a,long double b)
{
	if(a==(-numeric_limits<long double>::infinity()))
		return expFromLogValue(b);
	if(b==(-numeric_limits<long double>::infinity()))
		return expFromLogValue(a);

	long double m=max(a,b);
	return expFromLogValue(m+logl(expl(a-m)+expl(b-m)));
}

static bool computeAnchorEnvelope(kdAnchorAugNode*node,double*q,int dim,double& anchor_L,double& anchor_U,double*delta_u_out=NULL)
{
	if(node->anchor_W<=0 || !isfinite(node->anchor_W))
		return false;

	long double z2=0.0L;
	long double weighted_u_sum=0.0L;
	long double dot_min=0.0L;
	long double dot_max=0.0L;

	for(int d=0;d<dim;d++)
	{
		long double z=(long double)q[d]-(long double)node->anchor[d];
		long double y_lo=(long double)node->boundary[d][0]-(long double)node->anchor[d];
		long double y_hi=(long double)node->boundary[d][1]-(long double)node->anchor[d];
		long double prod_lo=z*y_lo;
		long double prod_hi=z*y_hi;

		z2+=z*z;
		weighted_u_sum+=2.0L*z*(long double)node->anchor_A[d];
		dot_min+=min(prod_lo,prod_hi);
		dot_max+=max(prod_lo,prod_hi);
	}

	long double W=(long double)node->anchor_W;
	long double logW=logl(W);
	long double mu=weighted_u_sum/W;
	long double u_min=2.0L*dot_min;
	long double u_max=2.0L*dot_max;
	long double z_norm=sqrtl(z2);
	long double ball_u_min=-2.0L*z_norm*(long double)node->anchor_radius;
	long double ball_u_max=2.0L*z_norm*(long double)node->anchor_radius;

	u_min=max(u_min,ball_u_min);
	u_max=min(u_max,ball_u_max);

	if(u_min>u_max)
		return false;

	if(delta_u_out!=NULL)
		*delta_u_out=(double)(u_max-u_min);

	if(mu<u_min)
		mu=u_min;
	if(mu>u_max)
		mu=u_max;

	anchor_L=expFromLogValue(logW-z2+mu);

	if(u_max-u_min<epsilon)
	{
		anchor_U=expFromLogValue(logW-z2+u_min);
			return isfinite(anchor_L) && isfinite(anchor_U);
	}

	long double denom=u_max-u_min;
	long double lambda_min=(u_max-mu)/denom;
	long double lambda_max=(mu-u_min)/denom;

	if(lambda_min<0.0L)
		lambda_min=0.0L;
	if(lambda_max<0.0L)
		lambda_max=0.0L;

	long double neg_inf=-numeric_limits<long double>::infinity();
	long double term_min=lambda_min>0.0L ? logl(lambda_min)+logW-z2+u_min : neg_inf;
	long double term_max=lambda_max>0.0L ? logl(lambda_max)+logW-z2+u_max : neg_inf;

	anchor_U=logSumExp2(term_min,term_max);

	if(anchor_U<anchor_L)
		anchor_U=anchor_L;

	return isfinite(anchor_L) && anchor_L>=0.0 && anchor_L<=node->sumE && anchor_U>=anchor_L;
}

static bool computeAdaptiveAnchorCandidate(kdAnchorAugNode*node,double*q,int dim,double& anchor_L,double& anchor_U)
{
	const double threshold=1.0;
	double delta_u;

	if(computeAnchorEnvelope(node,q,dim,anchor_L,anchor_U,&delta_u)==false)
		return false;

	double lb=ell_MBR(q,node->boundary,dim);
	double ub=u_MBR(q,node->boundary,dim);
	double delta_x=ub*ub-lb*lb;

	if(delta_x<epsilon)
		return false;

	return delta_u<threshold*delta_x;
}

double kdAnchorAugNode::LB(double*q,int dim,KDE_stat& stat)
{
	double base_L=kdLinearAugNode::LB(q,dim,stat);
	double anchor_L;
	double anchor_U;

	if(computeAnchorEnvelope(this,q,dim,anchor_L,anchor_U)==false)
		return base_L;

	return max(base_L,anchor_L);
}

double kdAnchorAugNode::UB(double*q,int dim,KDE_stat& stat)
{
	double base_U=kdLinearAugNode::UB(q,dim,stat);
	double anchor_L;
	double anchor_U;

	if(computeAnchorEnvelope(this,q,dim,anchor_L,anchor_U)==false)
		return base_U;

	return min(base_U,anchor_U);
}

kdAnchorAugNode*kdAnchorAugNode::createNode()
{
	return new kdAnchorAugNode();
}

double kdAnchorOnlyAugNode::LB(double*q,int dim,KDE_stat& stat)
{
	double anchor_L;
	double anchor_U;

	if(computeAnchorEnvelope(this,q,dim,anchor_L,anchor_U)==false)
		return kdLinearAugNode::LB(q,dim,stat);

	return anchor_L;
}

double kdAnchorOnlyAugNode::UB(double*q,int dim,KDE_stat& stat)
{
	double anchor_L;
	double anchor_U;

	if(computeAnchorEnvelope(this,q,dim,anchor_L,anchor_U)==false)
		return kdLinearAugNode::UB(q,dim,stat);

	return anchor_U;
}

kdAnchorOnlyAugNode*kdAnchorOnlyAugNode::createNode()
{
	return new kdAnchorOnlyAugNode();
}

double kdAdaptiveCombinedAnchorNode::LB(double*q,int dim,KDE_stat& stat)
{
	double base_L=kdLinearAugNode::LB(q,dim,stat);
	cache_q=q;
	cache_use_anchor=computeAdaptiveAnchorCandidate(this,q,dim,cache_anchor_L,cache_anchor_U);

	if(cache_use_anchor==false)
		return base_L;

	return max(base_L,cache_anchor_L);
}

double kdAdaptiveCombinedAnchorNode::UB(double*q,int dim,KDE_stat& stat)
{
	double base_U=kdLinearAugNode::UB(q,dim,stat);

	if(cache_q==q && cache_use_anchor==true)
		return min(base_U,cache_anchor_U);

	return base_U;
}

kdAdaptiveCombinedAnchorNode*kdAdaptiveCombinedAnchorNode::createNode()
{
	return new kdAdaptiveCombinedAnchorNode();
}

double kdAdaptiveSelectAnchorNode::LB(double*q,int dim,KDE_stat& stat)
{
	cache_q=q;
	cache_use_anchor=computeAdaptiveAnchorCandidate(this,q,dim,cache_anchor_L,cache_anchor_U);

	if(cache_use_anchor==true)
		return cache_anchor_L;

	return kdLinearAugNode::LB(q,dim,stat);
}

double kdAdaptiveSelectAnchorNode::UB(double*q,int dim,KDE_stat& stat)
{
	if(cache_q==q && cache_use_anchor==true)
		return cache_anchor_U;

	return kdLinearAugNode::UB(q,dim,stat);
}

kdAdaptiveSelectAnchorNode*kdAdaptiveSelectAnchorNode::createNode()
{
	return new kdAdaptiveSelectAnchorNode();
}

void kdAnchorAugNode::update_Aug(Node*node,Tree*t)
{
	this->updateAugInfo((kdLinearAugNode*)node,t);
	this->update_linearAugInfo((kdLinearAugNode*)node,t);
	this->update_anchorAugInfo((kdAnchorAugNode*)node,t);
}

void kdAnchorAugNode::update_anchorAugInfo(kdAnchorAugNode*node,Tree*t)
{
	int id;
	double y;
	double alpha;
	double y_norm_square;

	node->anchor=new double[t->dim];
	node->anchor_A=new double[t->dim];
	node->anchor_W=0;
	node->anchor_radius=0;

	for(int d=0;d<t->dim;d++)
	{
		node->anchor[d]=node->center[d];
		node->anchor_A[d]=0;
	}

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		y_norm_square=0;
		for(int d=0;d<t->dim;d++)
		{
			y=t->dataMatrix[id][d]-node->anchor[d];
			y_norm_square+=y*y;
		}
		if(y_norm_square>node->anchor_radius)
			node->anchor_radius=y_norm_square;

		alpha=exp(-y_norm_square);
		node->anchor_W+=alpha;

		for(int d=0;d<t->dim;d++)
		{
			y=t->dataMatrix[id][d]-node->anchor[d];
			node->anchor_A[d]+=alpha*y;
		}
	}
	node->anchor_radius=sqrt(node->anchor_radius);

	if((int)node->idList.size()<=t->leafCapacity)
		return;

	update_anchorAugInfo((kdAnchorAugNode*)node->childVector[0],t);
	update_anchorAugInfo((kdAnchorAugNode*)node->childVector[1],t);
}

//kd-tree
kdTree::kdTree(int dim,double**dataMatrix,int leafCapacity,KDE_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->leafCapacity=leafCapacity;
	this->stat=stat;
}

void kdTree::getNode_Boundary(kdNode*node)
{
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
		return;
	
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

void kdTree::build_kdTree(KDE_stat& stat)
{
	for(int i=0;i<stat.n;i++)
		rootNode->idList.push_back(i);

	getNode_Boundary((kdNode*)rootNode);
	KD_Tree_Recur((kdNode*)rootNode,0);
	initTree_sumE((kdNode*)rootNode);
}

void kdTree::initTree_sumE(kdNode*node)
{
	double sumE=0;
	int id;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		sumE+=1;
	}
	node->sumE=sumE;

	if((int)node->idList.size()>leafCapacity)
	{
		initTree_sumE((kdNode*)node->childVector[0]);
		initTree_sumE((kdNode*)node->childVector[1]);
	}
}

void kdTree::updateAugment(kdNode*node)
{
	node->update_Aug(node,this);
}
