#include "binaryTree.h"
#include "GBF_KC.h"

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
double binaryNode::LB(double*q,int dim,SVM_stat& stat)
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

double binaryNode::UB(double*q,int dim,SVM_stat& stat)
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
		return sum_alpha*exp(-stat.gammaValue*l2);

	exp_l2=exp(-stat.gammaValue*l2);
	exp_u2=exp(-stat.gammaValue*u2);
	//compute m and c
	m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*sum_alpha);
}

binaryNode*binaryNode::createNode()
{
	return new binaryNode();
}

void binaryNode::update_Aug(Node*node,Tree*t)
{
	//no code
}

void binaryNode::update_Aug(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	int id;
	double temp_norm_sq;
	double temp_radius;
	double max_radius=-inf;
	//initialization 
	sum_alpha=0;
	a_G=new double[dim];
	center=new double[dim];
	for(int d=0;d<dim;d++)
	{
		a_G[d]=0;
		center[d]=0;
	}
	S_G=0;

	//(1)update sum_alpha //(2)update a_G //(3)update S_G //(4)update center //(5)update radius
	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		sum_alpha+=alphaArray[id]; //(1)
		temp_norm_sq=0;
		for(int d=0;d<dim;d++)
		{
			a_G[d]+=alphaArray[id]*dataMatrix[id][d]; //(2)
			temp_norm_sq+=dataMatrix[id][d]*dataMatrix[id][d];
		}
		S_G+=alphaArray[id]*temp_norm_sq; //(3)
	}

	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sum_alpha; //(4)

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

void binaryNode::Aug_Incr(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat,int id/*,double& update_Radius*/)
{
	double square_norm=0;
	
	//(1)update sum_alpha //(2)update a_G //(3)update S_G //(4)update center //(5)update radius
	sum_alpha+=alphaArray[id]; //(1)
	for(int d=0;d<dim;d++)
	{
		a_G[d]+=alphaArray[id]*dataMatrix[id][d]; //(2)
		square_norm+=dataMatrix[id][d]*dataMatrix[id][d]; 
	}
	S_G+=alphaArray[id]*square_norm; //(3)
	
	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sum_alpha; //(4)
}

void binaryNode::Aug_Decr(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat,int id/*,double& update_Radius*/)
{
	double square_norm=0;
	//(1)update sum_alpha //(2)update a_G //(3)update S_G //(4)update center
	sum_alpha-=alphaArray[id]; //(1)
	for(int d=0;d<dim;d++)
	{
		a_G[d]-=alphaArray[id]*dataMatrix[id][d]; //(2)
		square_norm+=dataMatrix[id][d]*dataMatrix[id][d]; 
	}

	S_G-=alphaArray[id]*square_norm; //(3)
	
	if((int)idList.size()>0)
	{
		for(int d=0;d<dim;d++)
			center[d]=a_G[d]/sum_alpha; //(4)
	}
	else
	{
		for(int d=0;d<dim;d++)
			center[d]=0; //(4)
	}
}

void binaryNode::update_Radius(double**dataMatrix,double*alphaArray,int dim)
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

void half_division(vector<order_Entry>& o_EntryVector,binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	int halfSize=(int)ceil((double)o_EntryVector.size()/2.0);
	for(int i=0;i<(int)o_EntryVector.size();i++)
	{
		if(i<halfSize)
		{
			node1->idList.push_back(o_EntryVector[i].id);
			node1->Aug_Incr(dataMatrix,alphaArray,dim,stat,o_EntryVector[i].id);
		}
		else
		{
			node2->idList.push_back(o_EntryVector[i].id);
			node2->Aug_Incr(dataMatrix,alphaArray,dim,stat,o_EntryVector[i].id);
		}
	}

	node1->update_Radius(dataMatrix,alphaArray,dim);
	node2->update_Radius(dataMatrix,alphaArray,dim);

	o_EntryVector.clear();
}

void binaryNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector;
	//int halfSize;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		o_Entry.id=this->idList[i];
		o_Entry.dist=euclid_dist(dataMatrix[o_Entry.id],this->center,dim);
		o_EntryVector.push_back(o_Entry);
	}

	sort(o_EntryVector.begin(),o_EntryVector.end(),compare_order_Entry);

	half_division(o_EntryVector,node1,node2,dataMatrix,alphaArray,dim,stat);
}

void binaryNode::initNode(int dim)
{
	//(1)init sum_alpha //(2)init a_G //(3)init S_G //(4)init center //(5) radius
	sum_alpha=0; //(1)
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
	//(1)assign sum_alpha //(2)assign a_G //(3)assign S_G //(4)assign center //(5)assign radius
	sum_alpha=bNode_delta->sum_alpha; //(1)
	S_G=bNode_delta->S_G; //(3)
	//radius=bNode_delta->radius; //(5)
	for(int d=0;d<dim;d++)
	{
		a_G[d]=bNode_delta->a_G[d]; //(2)
		center[d]=bNode_delta->center[d]; //(4)
	}
}

void normNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		o_Entry.id=this->idList[i];
		o_Entry.norm=0;

		for(int d=0;d<dim;d++)
			o_Entry.norm+=dataMatrix[o_Entry.id][dim]*dataMatrix[o_Entry.id][dim];

		o_EntryVector.push_back(o_Entry);
	}

	sort(o_EntryVector.begin(),o_EntryVector.end(),compare_order_Entry_norm);
	half_division(o_EntryVector,node1,node2,dataMatrix,alphaArray,dim,stat);
}

normNode*normNode::createNode()
{
	return new normNode();
}

void two_far_ids(farNode*node,int& id1,int& id2,double**distMatrix,int dim)
{
	double max_dist=-inf;
	double curdist;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		for(int j=0;j<(int)node->idList.size();j++)
		{
			curdist=distMatrix[i][j];

			if(i!=j && curdist>max_dist)
			{
				max_dist=curdist;
				id1=node->idList[i];
				id2=node->idList[j];
			}
		}
	}

}

void createDistMatrix(double**& distMatrix,double**dataMatrix,int total_sv,int dim)
{
	double temp_sq_Euclid;

	distMatrix=new double*[total_sv];
	for(int i=0;i<total_sv;i++)
		distMatrix[i]=new double[total_sv];

	for(int i=0;i<total_sv;i++)
	{
		for(int j=0;j<total_sv;j++)
		{
			distMatrix[i][j]=0;
			for(int d=0;d<dim;d++)
			{
				temp_sq_Euclid=dataMatrix[i][d]-dataMatrix[j][d];
				distMatrix[i][j]+=temp_sq_Euclid*temp_sq_Euclid;
			}
		}
	}
}

void balance_saparation(double**distMatrix,int id1,int id2,vector<order_Entry>& o_EntryVector,binaryNode*node)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector1;
	vector<order_Entry> o_EntryVector2;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		o_Entry.id=node->idList[i];

		if((int)o_EntryVector1.size()>(int)(node->idList.size()/2.0))
		{
			o_EntryVector2.push_back(o_Entry);
			continue;
		}

		if((int)o_EntryVector2.size()>(int)(node->idList.size()/2.0))
		{
			o_EntryVector1.push_back(o_Entry);
			continue;
		}

		if(distMatrix[id1][node->idList[i]]<distMatrix[id2][node->idList[i]])
			o_EntryVector1.push_back(o_Entry);
		else
			o_EntryVector2.push_back(o_Entry);
	}

	o_EntryVector.insert(o_EntryVector.end(),o_EntryVector1.begin(),o_EntryVector1.end());
	o_EntryVector.insert(o_EntryVector.end(),o_EntryVector2.begin(),o_EntryVector2.end());
	o_EntryVector1.clear();
	o_EntryVector2.clear();
}

void farNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	vector<order_Entry> o_EntryVector;
	static double**distMatrix;
	static int isDistMatrixCreated=0;
	int id1,id2;

	if(isDistMatrixCreated==0)
	{
		createDistMatrix(distMatrix,dataMatrix,stat.total_sv,dim);
		isDistMatrixCreated++;
	}

	two_far_ids(this,id1,id2,distMatrix,dim);
	balance_saparation(distMatrix,id1,id2,o_EntryVector,this);
	half_division(o_EntryVector,node1,node2,dataMatrix,alphaArray,dim,stat);
}

farNode*farNode::createNode()
{
	return new farNode();
}

void ipNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	order_Entry o_Entry;
	vector<order_Entry> o_EntryVector;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		o_Entry.id=this->idList[i];
		o_Entry.ip=0;
		for(int d=0;d<dim;d++)
			o_Entry.ip+=dataMatrix[o_Entry.id][d]*this->center[d];
		o_EntryVector.push_back(o_Entry);
	}

	sort(o_EntryVector.begin(),o_EntryVector.end(),compare_order_Entry_ip);

	half_division(o_EntryVector,node1,node2,dataMatrix,alphaArray,dim,stat);
}

ipNode*ipNode::createNode()
{
	return new ipNode();
}

//ballNode
ballNode*ballNode::createNode()
{
	return new ballNode();
}

double ballNode::UB(double*q,int dim,SVM_stat& stat)
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

	//l_MBR and u_MBR
	//lb=ell_MBR(q,boundary,dim);
	//ub=u_MBR(q,boundary,dim);

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

void ballNode::updateBoundary(double**dataMatrix,SVM_stat& stat,int dim)
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

void ballNode::update_Aug(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	binaryNode::update_Aug(dataMatrix,alphaArray,dim,stat);
	updateBoundary(dataMatrix,stat,dim);
}

void find_furthest(double*pivot,vector<int>& idList,double**dataMatrix,double*alphaArray,int dim,int& id1,int& id2)
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

void ballNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	static int isRoot=1;
	int id1,id2;
	vector<order_Entry> o_EntryVector;

	if(isRoot==1)
		isRoot=0;
	else
		updateBoundary(dataMatrix,stat,dim);

	find_furthest(this->center,this->idList,dataMatrix,alphaArray,dim,id1,id2);
	balance_saparation_basic(dataMatrix,dim,id1,id2,o_EntryVector,this);

	half_division(o_EntryVector,node1,node2,dataMatrix,alphaArray,dim,stat);
}

ballNode_SOTA*ballNode_SOTA::createNode()
{
	return new ballNode_SOTA();
}

double ballNode_SOTA::LB(double*q,int dim,SVM_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sum_alpha*exp(-stat.gammaValue*ub*ub);

	return L;
}

double ballNode_SOTA::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sum_alpha*exp(-stat.gammaValue*lb*lb);

	return U;
}

oracleNode*oracleNode::createNode()
{
	return new oracleNode();
}

void oracleNode::GPA(binaryNode*node1,binaryNode*node2,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	int halfSize=(int)ceil((double)idList.size()/2.0);
	for(int i=0;i<(int)idList.size();i++)
	{
		if(i<halfSize)
		{
			node1->idList.push_back(idList[i]);
			node1->Aug_Incr(dataMatrix,alphaArray,dim,stat,idList[i]);
		}
		else
		{
			node2->idList.push_back(idList[i]);
			node2->Aug_Incr(dataMatrix,alphaArray,dim,stat,idList[i]);
		}
	}
}

double oracleNode::UB(double*q,int dim,SVM_stat& stat)
{
	double exp_l2,exp_u2;
	double l2,u2;
	double m,c;
	double dist;

	dist=euclid_dist(q,center,dim);

	if(innerRadius>dist || dist>outerRadius)
	{
		l2=max(innerRadius-dist,dist-outerRadius);
		l2=l2*l2;
	}
	else
		l2=0;

	u2=(outerRadius+dist)*(outerRadius+dist);

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

void oracleNode::updateOracleNode(double*ext,double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	int id;
	double temp_norm_sq;
	double temp_radius;
	double max_radius=-inf;
	double tempsq_ext_minus_p_Value;
	//initialization 
	sum_alpha=0;

	for(int d=0;d<dim;d++)
	{
		a_G[d]=0;
		center[d]=0;
	}
	S_G=0;

	innerRadius=inf;
	outerRadius=-inf;

	//(1)update sum_alpha //(2)update a_G //(3)update S_G //(4)update center //(5)update radius //(6) update innerRadius //(7) update outerRadius //(8) center
	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		sum_alpha+=alphaArray[id]; //(1)
		temp_norm_sq=0;
		for(int d=0;d<dim;d++)
		{
			a_G[d]+=alphaArray[id]*dataMatrix[id][d]; //(2)
			temp_norm_sq+=dataMatrix[id][d]*dataMatrix[id][d];
		}
		S_G+=alphaArray[id]*temp_norm_sq; //(3)
	}

	for(int d=0;d<dim;d++)
		center[d]=a_G[d]/sum_alpha; //(4)

	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		temp_radius=0;

		tempsq_ext_minus_p_Value=0;

		for(int d=0;d<dim;d++)
		{
			temp_radius+=(dataMatrix[id][d]-center[d])*(dataMatrix[id][d]-center[d]);
			tempsq_ext_minus_p_Value+=(ext[d]-dataMatrix[id][d])*(ext[d]-dataMatrix[id][d]);
		}

		if(temp_radius>max_radius)
			max_radius=temp_radius;

		if(innerRadius>tempsq_ext_minus_p_Value)
			innerRadius=tempsq_ext_minus_p_Value;
		if(outerRadius<tempsq_ext_minus_p_Value)
			outerRadius=tempsq_ext_minus_p_Value;
	}

	radius=sqrt(max_radius); //(5)
	innerRadius=sqrt(innerRadius); //(6)
	outerRadius=sqrt(outerRadius); //(7)
	//ring_center=ext; //(8)
	for(int d=0;d<dim;d++)
		center[d]=ext[d];
}

void oracleNode::updateOracle_LeafModel(double**dataMatrix,double*alphaArray,int dim,SVM_stat& stat)
{
	int elements=0;
	//reallocate the memory of the xspace and model in this leaf node
	for(int i=0;i<(int)idList.size();i++)
	{
		for(int d=0;d<dim;d++)
		{
			if(fabs(dataMatrix[idList[i]][d])>epsilon)
				elements++;
		}
	}

	elements+=(int)idList.size();

	x_space=(svm_node*)realloc(x_space,(elements)*sizeof(svm_node));
	model->sv_coef[0]=(double*)realloc(model->sv_coef[0],((int)idList.size())*sizeof(double));
	model->SV=(svm_node**)realloc(model->SV,((int)idList.size())*sizeof(svm_node*));

	//update model and x_space //(8) (9)
	createModel_inMemory(dataMatrix,alphaArray,dim,stat);
}

oracle_iNode*oracle_iNode::createNode()
{
	return new oracle_iNode();
}

double oracle_iNode::UB(double*q,int dim,SVM_stat& stat)
{
	return binaryNode::UB(q,dim,stat);
}

void binaryTree::build_BinaryTreeRecur(binaryNode*node)
{
	if((int)node->idList.size()>leafCapacity)
	{
		binaryNode*node1=node->createNode();
		binaryNode*node2=node->createNode();

		node1->initNode(dim);
		node2->initNode(dim);
		
		node->GPA(node1,node2,dataMatrix,alphaArray,dim,stat);

		build_BinaryTreeRecur(node1);
		build_BinaryTreeRecur(node2);

		node->childVector.push_back(node1);
		node->childVector.push_back(node2);
	}
	else //leafNode
	{
		node->update_Aug(dataMatrix,alphaArray,dim,stat);
		//with SVM model
		node->initModel_inMemory(dataMatrix,alphaArray,dim,stat);
		node->createModel_inMemory(dataMatrix,alphaArray,dim,stat);
	}
}

void binaryTree::build_BinaryTree()
{
	binaryNode*& bNode=(binaryNode*&) rootNode;
	for(int i=0;i<stat.total_sv;i++)
		bNode->idList.push_back(i);

	bNode->update_Aug(dataMatrix,alphaArray,dim,stat);

	build_BinaryTreeRecur((binaryNode*)rootNode);
}

binaryTree::binaryTree(int dim,double**dataMatrix,double*alphaArray,int leafCapacity,SVM_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->alphaArray=alphaArray;
	this->leafCapacity=leafCapacity;
	this->stat=stat;
}

void binaryTree::update_OracleRecur(double*ext,oracleNode*node)
{
	node->updateOracleNode(ext,dataMatrix,alphaArray,dim,stat);

	if((int)node->idList.size()>leafCapacity) //internal node
	{
		update_OracleRecur(ext,(oracleNode*)node->childVector[0]);
		update_OracleRecur(ext,(oracleNode*)node->childVector[1]);
	}
	else
		node->updateOracle_LeafModel(dataMatrix,alphaArray,dim,stat);
}

void binaryTree::update_Oracle(double*ext)
{
	oracleNode*& oNode=(oracleNode*&) rootNode;

	update_OracleRecur(ext,oNode);
}

void binaryTree::lookUp_BinaryTree_Recur(binaryNode*node)
{
	//(0)idList //(1)sum_alpha //(2)a_G //(3)S_G //(4)center //(5)radius
	cout<<"idList: "<<endl;
	for(int i=0;i<(int)node->idList.size();i++)
		cout<<node->idList[i]<<" "; //(0)
	cout<<endl;

	cout<<"sum_alpha: "<<node->sum_alpha<<endl; //(1)
	cout<<"a_G: ";
	for(int d=0;d<dim;d++)
		cout<<node->a_G[d]<<" "; //(2)
	cout<<endl;

	cout<<"S_G: "<<node->S_G<<endl; //(3)

	cout<<"center: ";
	for(int d=0;d<dim;d++)
		cout<<node->center[d]<<" "; //(4)
	cout<<endl;

	cout<<"radius: "<<node->radius<<endl; //(5)

	cout<<"__________________________________"<<endl;

	if((int)node->idList.size()>leafCapacity)
	{
		lookUp_BinaryTree_Recur((binaryNode*)node->childVector[0]);
		lookUp_BinaryTree_Recur((binaryNode*)node->childVector[1]);
	}
}

void binaryTree::lookUp_BinaryTree()
{
	lookUp_BinaryTree_Recur((binaryNode*)rootNode);
}

void ringTree::obtain_RingCenter()
{
	ring_center=new double[dim];

	double sum_alpha=0;
	for(int d=0;d<dim;d++)
		ring_center[d]=0;
	
	for(int i=0;i<stat.total_sv;i++)
	{
		sum_alpha+=alphaArray[i];
		for(int d=0;d<dim;d++)
			ring_center[d]+=alphaArray[i]*dataMatrix[i][d];
	}
	for(int d=0;d<dim;d++)
		ring_center[d]=ring_center[d]/sum_alpha;
}

void ringTree::build_RingTree()
{
	double**PP;
	double*reorder_alphaArray;
	vector<orderEntry> orderList;

	build_BinaryTree();

	obtain_RingCenter();
	init_reorder_mless(PP,reorder_alphaArray,stat,dim);
	oracle_Reorder(ring_center,dataMatrix,orderList,dim,stat);
	reordering_mless(dataMatrix,PP,alphaArray,reorder_alphaArray,orderList,stat,dim);
	orderList.clear();

	update_Oracle(ring_center);
}