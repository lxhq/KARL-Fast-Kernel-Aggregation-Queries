#include "m_tree.h"

mNode::mNode()
{
	parent=0;
}

double mNode::LB(double*q,int dim,SVM_stat& stat)
{
	double ub;
	double L;

	//u_tri
	ub=u_tri(q,O_r,dim,radius,temp_obt_dist);

	L=sum_alpha*exp(-stat.gammaValue*ub*ub);

	return L;
}

double mNode::UB(double*q,int dim,SVM_stat& stat)
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
	node->sum_alpha=0;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		for(int d=0;d<dim;d++)
			node->O_r[d]=node->O_r[d]+alphaArray[id]*dataMatrix[id][d];

		node->sum_alpha+=alphaArray[id];
	}
	for(int d=0;d<dim;d++)
		node->O_r[d]=node->O_r[d]/node->sum_alpha;

	//update radius
	node->radius=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		dist=euclid_dist(node->O_r,dataMatrix[id],dim);
		if(dist>node->radius)
			node->radius=dist;
	}

	//update model
	node->initModel_inMemory(dataMatrix,alphaArray,dim,stat);
	node->createModel_inMemory(dataMatrix,alphaArray,dim,stat);

	if(node->childVector.size()==0)
		return;

	for(int c=0;c<(int)node->childVector.size();c++)
		update_info((mNode*)node->childVector[c]);
}

//mAugNode
double mAugNode::LB(double*q,int dim,SVM_stat& stat)
{
	double ub;
	double L;

	ub=u_MBR(q,boundary,dim);

	L=sum_alpha*exp(-stat.gammaValue*ub*ub);

	return L;
}

double mAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double U;

	lb=ell_MBR(q,boundary,dim);

	U=sum_alpha*exp(-stat.gammaValue*lb*lb);

	return U;
}

//mLinearAugNode
double mLinearAugNode::LB(double*q,int dim,SVM_stat& stat)
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

double mLinearAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb,ub;
	double l2,u2;
	double exp_l2,exp_u2;
	double m,c;

	//l_tri and u_tri
	/*ub=u_tri(q,O_r,dim,radius,temp_obt_dist);
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
		return sum_alpha*exp(-stat.gammaValue*l2);

	exp_l2=exp(-stat.gammaValue*l2);
	exp_u2=exp(-stat.gammaValue*u2);
	//compute m and c
	m=(exp_u2-exp_l2)/(stat.gammaValue*(u2-l2));
	c=(u2*exp_l2-l2*exp_u2)/(u2-l2);

	return (m*gamma_sum+c*sum_alpha);
}

void mAugNode::updateBoundary(double**dataMatrix,SVM_stat& stat,int dim)
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

void mAugNode::update_Aug(Node*node,Tree*t)
{
	((mAugNode*)node)->updateBoundary(t->dataMatrix,t->stat,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

void mLinearAugNode::update_a_G(double**dataMatrix,double*alphaArray,int dim)
{
	a_G=new double[dim];
	int id;
	for(int d=0;d<dim;d++)
		a_G[d]=0;

	for(int i=0;i<(int)this->idList.size();i++)
	{
		id=idList[i];
		for(int d=0;d<dim;d++)
			a_G[d]+=alphaArray[id]*dataMatrix[id][d];
	}
}

void mLinearAugNode::update_S_G(double**dataMatrix,double*alphaArray,int dim)
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
		S_G+=alphaArray[id]*square_norm;
	}
}

void mLinearAugNode::update_Aug(Node*node,Tree*t)
{
	//update a_G and S_G
	((mLinearAugNode*)node)->update_a_G(t->dataMatrix,t->alphaArray,t->dim);
	((mLinearAugNode*)node)->update_S_G(t->dataMatrix,t->alphaArray,t->dim);
	((mLinearAugNode*)node)->updateBoundary(t->dataMatrix,t->stat,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

mTree::mTree(int dim,double**dataMatrix,double*alphaArray,int internalCapacity,int leafCapacity,SVM_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->alphaArray=alphaArray;
	this->internalCapacity=internalCapacity;
	this->leafCapacity=leafCapacity;
	this->stat=stat;

	O_p1=new double[dim];
	O_p2=new double[dim];
}

void mTree::insert(mNode*node,int id)
{
	if((int)node->childVector.size()==0) //leaf node
	{
		node->idList.push_back(id);
		if((int)node->idList.size()>leafCapacity)
			split(node,true);
	}
	else //non-leaf node
	{
		double expand_dist=inf;
		double inner_dist=inf;
		double dist;
		double diff;
		int child_id=-1;
		mNode* m_Node;
		for(int i=0;i<(int)node->childVector.size();i++)
		{
			m_Node=(mNode*)node->childVector[i];
			dist=euclid_dist(dataMatrix[id],m_Node->O_r,dim);
			if(dist<m_Node->radius)
			{
				inner_dist=min(inner_dist,dist);
				child_id=i;
			}
		}
		if(child_id==-1)
		{
			for(int i=0;i<(int)node->childVector.size();i++)
			{
				m_Node=(mNode*)node->childVector[i];
				diff=euclid_dist(dataMatrix[id],m_Node->O_r,dim)-m_Node->radius;
				if(diff<expand_dist)
				{
					expand_dist=min(expand_dist,diff);
					child_id=i;
				}
			}
			m_Node=(mNode*)node->childVector[child_id];
			m_Node->radius=euclid_dist(dataMatrix[id],m_Node->O_r,dim);
		}

		insert((mNode*)node->childVector[child_id],id);
	}
}

void mTree::split(mNode*node,bool isLeaf)
{
	//mNode*node_dash=new mNode();
	mNode*node_dash=node->createNode();
	node_dash->O_r=new double[dim];

	promote(node,isLeaf);
	partition(node,node_dash,isLeaf);

	//node is the root
	if(node->parent==0)
	{
		mNode*newRoot=node->createNode();
		//mNode*newRoot=new mNode();
		newRoot->O_r=new double[dim];
		node->parent=newRoot;
		node_dash->parent=newRoot;
		newRoot->childVector.push_back(node);
		newRoot->childVector.push_back(node_dash);
		rootNode=newRoot;
	}
	else
	{
		mNode*parentNode=node->parent;
		node_dash->parent=parentNode;
		parentNode->childVector.push_back(node_dash);
		if((int)parentNode->childVector.size()>internalCapacity)
			split(parentNode,false);
	}
}

void mTree::promote(mNode*node,bool isLeaf)
{
	int size;
	if(isLeaf==true)
		size=(int)node->idList.size();
	else
		size=(int)node->childVector.size();

	distMatrix=new double*[size];
	for(int i=0;i<size;i++)
		distMatrix[i]=new double[size];

	if(isLeaf==true)
	{
		for(int i=0;i<(int)node->idList.size();i++)
			for(int j=0;j<(int)node->idList.size();j++)
				distMatrix[i][j]=euclid_dist(dataMatrix[node->idList[i]],dataMatrix[node->idList[j]],dim);
	}
	else
	{
		for(int i=0;i<(int)node->childVector.size();i++)
			for(int j=0;j<(int)node->childVector.size();j++)
				distMatrix[i][j]=euclid_dist(((mNode*)node->childVector[i])->O_r,((mNode*)node->childVector[j])->O_r,dim);
	}

	int best_i,best_j;
	double radius_i,radius_j;
	double min_radius=inf;
	double max_radius;
	for(int i=0;i<size;i++)
	{
		for(int j=i+1;j<size;j++)
		{
			radius_i=0;
			radius_j=0;

			for(int k=0;k<size;k++)
			{
				if(distMatrix[i][k]<=distMatrix[j][k])
					radius_i=max(radius_i,distMatrix[i][k]);
				else
					radius_j=max(radius_j,distMatrix[j][k]);

				max_radius=max(radius_i,radius_j);
			}

			if(max_radius<min_radius)
			{
				min_radius=max_radius;
				best_i=i;
				best_j=j;
			}

		}
	}

	if(isLeaf==true)
	{
		for(int d=0;d<dim;d++)
		{
			O_p1[d]=dataMatrix[node->idList[best_i]][d];
			O_p2[d]=dataMatrix[node->idList[best_j]][d];
		}
	}
	else
	{
		for(int d=0;d<dim;d++)
		{
			O_p1[d]=((mNode*)node->childVector[best_i])->O_r[d];
			O_p2[d]=((mNode*)node->childVector[best_j])->O_r[d];
		}
	}

	for(int i=0;i<size;i++)
		delete[] distMatrix[i];
	delete distMatrix;
}

void mTree::partition(mNode*node,mNode*node_dash,bool isLeaf)
{
	double radius_one;
	double dist_one;
	double dist_two;
	vector<int> idNodeList;
	vector<mNode*> nodeChildVector; 

	radius_one=0;
	node_dash->radius=0;
	if(isLeaf==true)
	{
		int id;
		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			dist_one=euclid_dist(O_p1,dataMatrix[id],dim);
			dist_two=euclid_dist(O_p2,dataMatrix[id],dim);
			if(dist_one<=dist_two && (int)idNodeList.size()<=(int)floor((double)leafCapacity/2.0))
			{
				idNodeList.push_back(id);
				radius_one=max(radius_one,dist_one);
			}
			else
			{
				node_dash->idList.push_back(id);
				node_dash->radius=max(node_dash->radius,dist_two);
			}
		}

		node->idList.clear();
		for(int i=0;i<(int)idNodeList.size();i++)
			node->idList.push_back(idNodeList[i]);
	}
	else //non-leaf node
	{
		mNode*m_Node;
		for(int c=0;c<(int)node->childVector.size();c++)
		{
			m_Node=(mNode*)node->childVector[c];
			dist_one=euclid_dist(O_p1,m_Node->O_r,dim);
			dist_two=euclid_dist(O_p2,m_Node->O_r,dim);
			if(dist_one<=dist_two && (int)nodeChildVector.size()<=(int)floor((double)internalCapacity/2.0))
			{
				nodeChildVector.push_back((mNode*)node->childVector[c]);
				radius_one=max(radius_one,dist_one+m_Node->radius);
			}
			else
			{
				node_dash->childVector.push_back(node->childVector[c]);
				node_dash->radius=max(node_dash->radius,dist_two+m_Node->radius);
			}
		}
		node->childVector.clear();
		for(int c=0;c<(int)nodeChildVector.size();c++)
			node->childVector.push_back(nodeChildVector[c]);
	}

	//Update Radius in node and O_r in two child nodes
	node->radius=radius_one;
	for(int d=0;d<dim;d++)
	{
		node->O_r[d]=O_p1[d];
		node_dash->O_r[d]=O_p2[d];
	}
}

void mTree::build_m_tree()
{
	//initialization of O_r for root
	((mNode*)rootNode)->O_r=new double[dim];

	for(int i=0;i<stat.total_sv;i++)
		insert((mNode*)rootNode,i);

	//Augment info
	update_Aug((mNode*)rootNode);
	update_rootInfo();
}

void mTree::lookUp_m_tree()
{
	if(rootNode->childVector.size()==0)//rootNode is the child node
	{
		cout<<"rootNode idList: ";
		for(int i=0;i<(int)rootNode->idList.size();i++)
			cout<<rootNode->idList[i]<<" ";
		cout<<endl;
	}

	for(int c=0;c<(int)rootNode->childVector.size();c++)
		lookUp_m_tree_Recur((mNode*)rootNode->childVector[c]);
}

void mTree::lookUp_m_tree_Recur(mNode*node)
{
	cout<<"radius: "<<node->radius<<endl;
	cout<<"O_r"<<endl;
	for(int d=0;d<dim;d++)
		cout<<node->O_r[d]<<" ";

	cout<<"idList: ";
	for(int i=0;i<(int)node->idList.size();i++)
		cout<<node->idList[i]<<" ";
	cout<<endl;

	cout<<"sum_alpha: "<<node->sum_alpha<<endl;

	cout<<endl;

	for(int c=0;c<(int)node->childVector.size();c++)
		lookUp_m_tree_Recur((mNode*)node->childVector[c]);
}

void mTree::update_Aug(mNode*node)
{
	int id;
	if(node->childVector.size()==0)
	{
		node->sum_alpha=0;
		for(int i=0;i<(int)node->idList.size();i++)
		{
			id=node->idList[i];
			node->sum_alpha+=this->alphaArray[id];
		}

		node->initModel_inMemory(dataMatrix,alphaArray,dim,stat);
		node->createModel_inMemory(dataMatrix,alphaArray,dim,stat);
	}
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug((mNode*)node->childVector[c]);

	for(int c=0;c<(int)node->childVector.size();c++)
		for(int i=0;i<(int)node->childVector[c]->idList.size();i++)
			node->idList.push_back(node->childVector[c]->idList[i]);

	node->sum_alpha=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		node->sum_alpha+=this->alphaArray[id];
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
	for(int i=0;i<stat.total_sv;i++)
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
	double build_Time;
	treeFile.open(treeFileName);

	if(treeFile.is_open()==false)
	{
		cout<<"Tree has not been built!"<<endl;

		#ifdef C_PLUSPLUS11_CLOCK
			auto start_b_index_s=chrono::high_resolution_clock::now();
		#else
			clock_t start_b_index_s=clock();
		#endif

		build_BL_m_tree();

		#ifdef C_PLUSPLUS11_CLOCK
			auto end_b_index_s=chrono::high_resolution_clock::now();
			build_Time=(chrono::duration_cast<chrono::nanoseconds>(end_b_index_s-start_b_index_s).count())/1000000000.0;
		#else
			clock_t end_b_index_s=clock();
			build_Time=((double)(end_b_index_s-start_b_index_s))/CLOCKS_PER_SEC;
		#endif

		cout<<"Build m-tree Time is: "<<build_Time<<"sec"<<endl;
		
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