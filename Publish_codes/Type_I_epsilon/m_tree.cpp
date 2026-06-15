#include "m_tree.h"

double mNode::LB(double*q,int dim,KDE_stat& stat)
{
	double ub;
	double L;

	//u_tri
	ub=u_tri(q,O_r,dim,radius,temp_obt_dist);

	L=sumE*exp(-ub*ub);

	return L;
}

double mNode::UB(double*q,int dim,KDE_stat& stat)
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

//mAugNode
double mAugNode::LB(double*q,int dim,KDE_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sumE*exp(-ub*ub);

	return L;
}

double mAugNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sumE*exp(-lb*lb);

	return U;
}

/*void mAugNode::updateBoundary(double**dataMatrix,KDE_stat& stat,int dim)
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
}*/

void mAugNode::update_Aug(Node*node,Tree*t)
{
	//((mAugNode*)node)->updateBoundary(t->dataMatrix,t->stat,t->dim);
	updateBoundary(t->dataMatrix,node->idList,((mAugNode*)node)->boundary,t->stat,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

//mLinearAugNode
double mLinearAugNode::LB(double*q,int dim,KDE_stat& stat)
{
	double t_star;
	double ip=0;
	for(int d=0;d<dim;d++)
		ip=ip+q[d]*a_G[d];

	//Compute t^*
	gamma_sum=sumE*stat.qSquareNorm-2*ip+S_G;
	t_star=(gamma_sum/sumE);

	return sumE*exp(-t_star);
}

double mLinearAugNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//l_tri and u_tri
	ub=u_tri(q,O_r,dim,radius,temp_obt_dist);
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

void mLinearAugNode::update_a_G(double**dataMatrix,int dim)
{
	a_G=new double[dim];
	int id;
	for(int d=0;d<dim;d++)
		a_G[d]=0;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		id=idList[i];
		for(int d=0;d<dim;d++)
			a_G[d]+=dataMatrix[id][d];
	}
}

void mLinearAugNode::update_S_G(double**dataMatrix,int dim)
{
	int id;
	S_G=0;
	double square_norm;
	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		square_norm=0;
		for(int d=0;d<dim;d++)
			square_norm+=dataMatrix[id][d]*dataMatrix[id][d];
		S_G+=square_norm;
	}
}

void mLinearAugNode::update_Aug(Node*node,Tree*t)
{
	//update a_G and S_G
	((mLinearAugNode*)node)->update_a_G(t->dataMatrix,t->dim);
	((mLinearAugNode*)node)->update_S_G(t->dataMatrix,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

double mLinearAug_RectNode::UB(double*q,int dim,KDE_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

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

void mLinearAug_RectNode::update_Aug(Node *node, Tree *t)
{
	//update a_G, S_G and boundary
	((mLinearAug_RectNode*)node)->update_a_G(t->dataMatrix,t->dim);
	((mLinearAug_RectNode*)node)->update_S_G(t->dataMatrix,t->dim);
	updateBoundary(t->dataMatrix,node->idList,((mLinearAug_RectNode*)node)->boundary,t->stat,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

//mTree
mTree::mTree(int dim,double**dataMatrix,int internalCapacity,int leafCapacity,KDE_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->internalCapacity=internalCapacity;
	this->leafCapacity=leafCapacity;
	this->stat=stat;

	O_p1=new double[dim];
	O_p2=new double[dim];
}

void mTree::update_Aug(mNode*node)
{
	int id;
	if(node->childVector.size()==0)
	{
		node->sumE=0;
		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			node->sumE+=1;
		}
	}
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug((mNode*)node->childVector[c]);

	for(int c=0;c<(int)node->childVector.size();c++)
		for(int i=0;i<(int)node->childVector[c]->idList.size();i++)
			node->idList.push_back(node->childVector[c]->idList[i]);

	node->sumE=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		node->sumE+=1;
	}
}

void mTree::update_rootInfo()
{
	mNode*& root=(mNode*&)this->rootNode;
	if(root->childVector.size()==0)
		return;
	//update center
	for(int d=0;d<dim;d++)
		root->O_r[d]=((mNode*)root->childVector[0])->O_r[d];

	//update radius
	root->radius=-inf;
	double temp_Radius;
	for(int c=0;c<(int)root->childVector.size();c++)
	{
		temp_Radius=euclid_dist(root->O_r,((mNode*)root->childVector[c])->O_r,dim)+((mNode*)root->childVector[c])->radius;
		if(temp_Radius>root->radius)
			root->radius=temp_Radius;
	}
}

void mTree::updateAugment(mNode*node)
{
	node->update_Aug(node,this);
}

int which_bestPivot(int id,double**dataMatrix,vector<int>& pivot_idList,int dim,vector<bool>& _is_pivot)
{
	double best_dist=inf;
	int best_pivot=-1;
	double dist;

	for(int p=0;p<(int)pivot_idList.size();p++)
	{
		if(_is_pivot[p]==false)
			continue;

		dist=euclid_dist(dataMatrix[pivot_idList[p]],dataMatrix[id],dim);

		if(dist<best_dist)
		{
			best_dist=dist;
			best_pivot=p;
		}
	}
	return best_pivot;
}

void re_init_List(vector< vector<int> >& _idList,vector<int>& pivot_idList,vector<bool>& _is_pivot)
{
	for(int p=0;p<(int)pivot_idList.size();p++)
	{
		_is_pivot[p]=true;
		_idList[p].clear();
	}

	pivot_idList.clear();
}

void find_furthest(double**dataMatrix,mNode*node,int dim,int& report_id1,int& report_id2)
{
	double dist_max=-inf;
	double dist;
	int id1;
	int id2;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id1=node->idList[i];
		for(int j=i+1;j<(int)node->idList.size();j++)
		{
			id2=node->idList[j];
			dist=euclid_dist(dataMatrix[id1],dataMatrix[id2],dim);
			if(dist>dist_max)
			{
				dist_max=dist;
				report_id1=id1;
				report_id2=id2;
			}
		}
	}
}

void binary_partition(mNode*node,vector<int>& pivot_idList,vector< vector<int> >& _idList,double**dataMatrix,int dim)
{
	int id1;
	int id2;
	int id;
	double dist1;
	double dist2;
	int div_size=(int)(node->idList.size()/2.0);

	find_furthest(dataMatrix,node,dim,id1,id2);
	pivot_idList.push_back(id1);
	pivot_idList.push_back(id2);

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(id==id1)
		{
			_idList[0].push_back(id);
			continue;
		}
		if(id==id2)
		{
			_idList[1].push_back(id);
			continue;
		}

		if((int)_idList[0].size()>div_size)
		{
			_idList[1].push_back(id);
			continue;
		}
		if((int)_idList[1].size()>div_size)
		{
			_idList[0].push_back(id);
			continue;
		}

		dist1=euclid_dist(dataMatrix[id1],dataMatrix[id],dim);
		dist2=euclid_dist(dataMatrix[id2],dataMatrix[id],dim);

		if(dist1<dist2)
			_idList[0].push_back(id);
		else
			_idList[1].push_back(id);
	}
}

void mTree::build_BL_m_tree_Recur(mNode*node)
{
	vector< vector<int> > _idList;
	vector<int> temp_idList;
	vector<int> pivot_idList;
	vector<bool> _is_pivot;
	int best_pivot;
	int id;
	int true_Counter=0;
	int counter=0;
	int loop_Counter=0;

	int m=max(internalCapacity/2,(int)(node->idList.size()/(2.0*(double)internalCapacity)));

	if((int)node->idList.size()<=leafCapacity)
		return;

	for(int i=0;i<internalCapacity;i++)
	{
		_is_pivot.push_back(true);
		_idList.push_back(temp_idList);
	}
	
	while(true_Counter<=1)
	{
		true_Counter=0;

		re_init_List(_idList,pivot_idList,_is_pivot);

		//There are too many loops, this group seems to be compact
		if(loop_Counter>=5)
		{
			//use binary-partition method to avoid many loops
			binary_partition(node,pivot_idList,_idList,dataMatrix,dim);

			true_Counter=2;
			break;
		}
		else
		{
			re_init_List(_idList,pivot_idList,_is_pivot);
			sample(node,pivot_idList,internalCapacity);
		}

		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			
			best_pivot=which_bestPivot(id,dataMatrix,pivot_idList,dim,_is_pivot);

			_idList[best_pivot].push_back(id);
		}

		//delete sub-groups which contain too few objects
		for(int p=0;p<(int)pivot_idList.size();p++)
		{
			if((int)_idList[p].size()<m)
			{
				_is_pivot[p]=false;
				for(int i=0;i<(int)_idList[p].size();i++)
				{
					best_pivot=which_bestPivot(_idList[p][i],dataMatrix,pivot_idList,dim,_is_pivot);
					_idList[best_pivot].push_back(_idList[p][i]);
				}
				_idList[p].clear();
			}
		}

		for(int p=0;p<(int)pivot_idList.size();p++)
		{
			if(_is_pivot[p]==true)
				true_Counter++;
		}		

		loop_Counter++;
	}

	for(int c=0;c<true_Counter;c++)
		node->childVector.push_back((Node*)node->createNode());

	counter=0;
	for(int p=0;p<(int)pivot_idList.size();p++)
	{
		if(_is_pivot[p]==true)
		{
			for(int i=0;i<(int)_idList[p].size();i++)
				node->childVector[counter]->idList.push_back(_idList[p][i]);
			counter++;
		}
	}
	
	for(int c=0;c<(int)node->childVector.size();c++)
		build_BL_m_tree_Recur((mNode*)node->childVector[c]);
}

void mTree::build_BL_m_tree()
{
	for(int i=0;i<stat.n;i++)
		rootNode->idList.push_back(i);

	build_BL_m_tree_Recur((mNode*)rootNode);
	update_BL_Augment((mNode*)rootNode);
}

void mTree::update_BL_Augment(mNode*node)
{
	update_info(node);
	node->update_Aug(node,this);
}

void mTree::sample(mNode*node,vector<int>& sample_idList,int sampleNum)
{
	//Initialization in the sample function
	static int accessFunc=0;
	if(accessFunc==0)
	{
		//srand(1);
		srand(time(0));
		accessFunc++;
	}

	//Sampling
	vector<int> indexVector;
	int size=node->idList.size();
	
	for(int i=0;i<size;i++)
		indexVector.push_back(i);

	random_shuffle(indexVector.begin(),indexVector.end());

	for(int s=0;s<sampleNum;s++)
		sample_idList.push_back(node->idList[indexVector[s]]);
}

//This function is similar with update_Aug in mTree class which is used in bulk-loading
void mTree::update_info(mNode*node)
{
	//no need to update parent in bulk-loading algorithm
	int id;
	double dist;

	//update O_r (O_r is used as the center of the group)
	node->O_r=new double[dim];
	for(int d=0;d<dim;d++)
		node->O_r[d]=0;
	node->sumE=0;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		for(int d=0;d<dim;d++)
			node->O_r[d]=node->O_r[d]+dataMatrix[id][d];

		node->sumE+=1;
	}
	for(int d=0;d<dim;d++)
		node->O_r[d]=node->O_r[d]/node->sumE;

	//update radius
	node->radius=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		dist=euclid_dist(node->O_r,dataMatrix[id],dim);
		if(dist>node->radius)
			node->radius=dist;
	}

	if(node->childVector.size()==0)
		return;

	for(int c=0;c<(int)node->childVector.size();c++)
		update_info((mNode*)node->childVector[c]);
}

void save_Tree_Recur(fstream& treeFile,mNode*node)
{
	treeFile<<(int)node->idList.size()<<endl;
	for(int i=0;i<(int)node->idList.size();i++)
		treeFile<<node->idList[i]<<" ";
	treeFile<<endl;
	treeFile<<(int)node->childVector.size()<<endl;

	for(int c=0;c<(int)node->childVector.size();c++)
		save_Tree_Recur(treeFile,(mNode*)node->childVector[c]);
}

void load_Tree_Recur(fstream& treeFile,mNode*node)
{
	int idList_size;
	int id;
	int childSize;
	treeFile>>idList_size;

	for(int i=0;i<idList_size;i++)
	{
		treeFile>>id;
		node->idList.push_back(id);
	}
	treeFile>>childSize;

	for(int c=0;c<childSize;c++)
		node->childVector.push_back((mNode*)node->createNode());

	for(int c=0;c<childSize;c++)
		load_Tree_Recur(treeFile,(mNode*)node->childVector[c]);
}

void mTree::save_Tree(char*treeFileName)
{
	fstream treeFile;
	treeFile.open(treeFileName,ios::in | ios::out | ios::trunc);
	if(treeFile.is_open()==false)
	{
		cout<<"Cannot Save Tree!"<<endl;
		exit(1);
	}

	save_Tree_Recur(treeFile,(mNode*)rootNode);

	treeFile.close();
	cout<<"Finish saving the tree!"<<endl;
}

void mTree::load_Tree(char*treeFileName)
{
	fstream treeFile;
	treeFile.open(treeFileName);

	if(treeFile.is_open()==false)
	{
		cout<<"Tree has not been built!"<<endl;
		build_BL_m_tree();
		save_Tree(treeFileName);
		return;
	}
	else
	{
		load_Tree_Recur(treeFile,(mNode*)rootNode);
		update_BL_Augment((mNode*)rootNode);
	}

	treeFile.close();
}

void updateBoundary(double**dataMatrix,vector<int>& idList,double**& boundary,KDE_stat& stat,int dim)
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