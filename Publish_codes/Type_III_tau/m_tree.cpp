#include "m_tree.h"

mNode::mNode()
{
	parent=0;
}

double mNode::LB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double ub;
	double L_pos;
	double U_neg;

	if(sum_alpha_pos<epsilon)
		ub=0;
	else
		ub=u_tri(q,O_r_pos,dim,radius_pos,temp_obt_dist_pos); //G_{+}

	if(sum_alpha_neg<epsilon)
		lb=0;
	else
		lb=ell_tri(q,O_r_neg,dim,radius_neg,temp_obt_dist_neg); //G_{-}

	L_pos=sum_alpha_pos*exp(-stat.gammaValue*ub*ub); //LB_{G_{+}}
	U_neg=sum_alpha_neg*exp(-stat.gammaValue*lb*lb); //UB_{G_{-}}

	return (L_pos-U_neg);
}

double mNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double ub;
	double L_neg;
	double U_pos;

	//G_{+}
	if(sum_alpha_pos<epsilon)
		lb=0;
	else
	{
		if(temp_obt_dist_pos<radius_pos)
			lb=0;
		else
			lb=temp_obt_dist_pos-radius_pos;
	}
	
	//G_{-}
	if(sum_alpha_neg<epsilon)
		ub=0;
	else
		ub=temp_obt_dist_neg+radius_neg;

	L_neg=sum_alpha_neg*exp(-stat.gammaValue*ub*ub); //LB_{G_{-}}
	U_pos=sum_alpha_pos*exp(-stat.gammaValue*lb*lb); //UB_{G_{+}}

	return (U_pos-L_neg);
}

//mAugNode
double mAugNode::LB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double ub;
	double L_pos;
	double U_neg;

	if(sum_alpha_pos<epsilon)
		ub=0;
	else
		ub=u_MBR(q,boundary_pos,dim); //G_{+}	

	if(sum_alpha_neg<epsilon)
		lb=0;
	else
		lb=ell_MBR(q,boundary_neg,dim); //G_{-}

	L_pos=sum_alpha_pos*exp(-stat.gammaValue*ub*ub); //LB_{G_{+}}
	U_neg=sum_alpha_neg*exp(-stat.gammaValue*lb*lb); //UB_{G_{-}}

	return (L_pos-U_neg);
}

double mAugNode::UB(double*q,int dim,SVM_stat& stat)
{
	double lb;
	double ub;
	double L_neg;
	double U_pos;

	if(sum_alpha_neg<epsilon)
		ub=0;
	else
		ub=u_MBR(q,boundary_neg,dim); //G_{-}

	if(sum_alpha_pos<epsilon)
		lb=0;
	else
		lb=ell_MBR(q,boundary_pos,dim); //G_{+}

	L_neg=sum_alpha_neg*exp(-stat.gammaValue*ub*ub); //LB_{G_{-}}
	U_pos=sum_alpha_pos*exp(-stat.gammaValue*lb*lb); //UB_{G_{+}}

	return (U_pos-L_neg);
}

void mAugNode::updateBoundary(double**dataMatrix,double*outputArray,SVM_stat& stat,int dim)
{
	int id;

	boundary_pos=new double*[dim];
	boundary_neg=new double*[dim];
	for(int d=0;d<dim;d++)
	{
		boundary_pos[d]=new double[2];
		boundary_neg[d]=new double[2];
	}

	//initialization of boundary
	for(int d=0;d<dim;d++)
	{
		boundary_pos[d][0]=inf;
		boundary_pos[d][1]=-inf;

		boundary_neg[d][0]=inf;
		boundary_neg[d][1]=-inf;
	}

	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];

		if(outputArray[id]>=0)
		{
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[id][d]<boundary_pos[d][0])
					boundary_pos[d][0]=dataMatrix[id][d];
				if(dataMatrix[id][d]>boundary_pos[d][1])
					boundary_pos[d][1]=dataMatrix[id][d];
			}
		}
		else
		{
			for(int d=0;d<dim;d++)
			{
				if(dataMatrix[id][d]<boundary_neg[d][0])
					boundary_neg[d][0]=dataMatrix[id][d];
				if(dataMatrix[id][d]>boundary_neg[d][1])
					boundary_neg[d][1]=dataMatrix[id][d];
			}
		}
	}
}

void mAugNode::update_Aug(Node*node,Tree*t)
{
	((mAugNode*)node)->updateBoundary(t->dataMatrix,t->outputArray,t->stat,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

//mLinearAugNode
double mLinearAugNode::LB(double*q,int dim,SVM_stat& stat)
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
		ub=u_tri(q,O_r_neg,dim,radius_neg,temp_obt_dist_neg);
		if(temp_obt_dist_neg<radius_neg)
			lb=0;
		else
			lb=temp_obt_dist_neg-radius_neg;

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

double mLinearAugNode::UB(double*q,int dim,SVM_stat& stat)
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
		ub=u_tri(q,O_r_pos,dim,radius_pos,temp_obt_dist_pos);
		if(temp_obt_dist_pos<radius_pos)
			lb=0;
		else
			lb=temp_obt_dist_pos-radius_pos;

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

void mLinearAugNode::update_a_G(double**dataMatrix,double*outputArray,int dim)
{
	a_G_pos=new double[dim];
	a_G_neg=new double[dim];
	int id;

	for(int d=0;d<dim;d++)
	{
		a_G_pos[d]=0;
		a_G_neg[d]=0;
	}

	for(int i=0;i<(int)this->idList.size();i++)
	{
		id=idList[i];

		if(outputArray[id]>=0)
		{
			for(int d=0;d<dim;d++)
				a_G_pos[d]+=outputArray[id]*dataMatrix[id][d];
		}
		else
		{
			for(int d=0;d<dim;d++)
				a_G_neg[d]+=fabs(outputArray[id])*dataMatrix[id][d];
		}	
	}
}

void mLinearAugNode::update_S_G(double**dataMatrix,double*outputArray,int dim)
{
	int id;
	S_G_pos=0;
	S_G_neg=0;

	double square_norm;
	for(int i=0;i<(int)idList.size();i++)
	{
		id=idList[i];
		square_norm=0;
		for(int d=0;d<dim;d++)
			square_norm+=dataMatrix[id][d]*dataMatrix[id][d];

		if(outputArray[id]>=0)
			S_G_pos+=outputArray[id]*square_norm;
		else
			S_G_neg+=fabs(outputArray[id])*square_norm;
	}
}

void mLinearAugNode::update_Aug(Node*node,Tree*t)
{
	//update a_G and S_G
	((mLinearAugNode*)node)->update_a_G(t->dataMatrix,t->outputArray,t->dim);
	((mLinearAugNode*)node)->update_S_G(t->dataMatrix,t->outputArray,t->dim);
	for(int c=0;c<(int)node->childVector.size();c++)
		update_Aug(node->childVector[c],t);
}

//Tree part
mTree::mTree(int dim,double**dataMatrix,double*outputArray,int internalCapacity,int leafCapacity,SVM_stat& stat)
{
	this->dim=dim;
	this->dataMatrix=dataMatrix;
	this->outputArray=outputArray;
	this->internalCapacity=internalCapacity;
	this->leafCapacity=leafCapacity;
	this->stat=stat;

	O_p1=new double[dim];
	O_p2=new double[dim];
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

//This function is similar with update_Aug in mTree class which is used in bulk-loading
void mTree::update_info(mNode*node)
{
	//no need to update parent in bulk-loading algorithm
	int id;
	double dist;

	//update O_r (O_r is used as the center of the group)
	node->O_r_pos=new double[dim];
	node->O_r_neg=new double[dim];
	for(int d=0;d<dim;d++)
	{
		node->O_r_pos[d]=0;
		node->O_r_neg[d]=0;
	}
	node->sum_alpha_pos=0;
	node->sum_alpha_neg=0;

	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		for(int d=0;d<dim;d++)
		{
			if(outputArray[id]>=0)
				node->O_r_pos[d]=node->O_r_pos[d]+outputArray[id]*dataMatrix[id][d];
			else
				node->O_r_neg[d]=node->O_r_neg[d]+fabs(outputArray[id])*dataMatrix[id][d];
		}

		if(outputArray[id]>=0)
			node->sum_alpha_pos+=outputArray[id];
		else
			node->sum_alpha_neg+=fabs(outputArray[id]);
	}

	for(int d=0;d<dim;d++)
	{
		node->O_r_pos[d]=node->O_r_pos[d]/node->sum_alpha_pos;
		node->O_r_neg[d]=node->O_r_neg[d]/node->sum_alpha_neg;
	}

	//update radius
	node->radius_pos=0;
	node->radius_neg=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		id=node->idList[i];
		if(outputArray[id]>=0)
		{
			dist=euclid_dist(node->O_r_pos,dataMatrix[id],dim);
			if(dist>node->radius_pos)
				node->radius_pos=dist;
		}
		else
		{
			dist=euclid_dist(node->O_r_neg,dataMatrix[id],dim);
			if(dist>node->radius_neg)
				node->radius_neg=dist;
		}
	}

	//update model
	node->initModel_inMemory(dataMatrix,outputArray,dim,stat);
	node->createModel_inMemory(dataMatrix,outputArray,dim,stat);

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