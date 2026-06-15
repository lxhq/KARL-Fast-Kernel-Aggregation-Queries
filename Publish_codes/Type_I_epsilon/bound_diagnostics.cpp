#include "init_KAQ.h"
#include "kd_tree.h"

#include <iomanip>

static double sqNorm(double*q,int dim)
{
	double value=0;
	for(int d=0;d<dim;d++)
		value+=q[d]*q[d];

	return value;
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

struct AnchorDiag
{
	bool valid;
	double lb;
	double ub;
	double delta_u;
	double u_min;
	double u_max;
	double z_norm;
	double radius;
};

static AnchorDiag computeAnchorDiag(kdAnchorAugNode*node,double*q,int dim)
{
	AnchorDiag diag={};
	diag.valid=false;
	diag.lb=0;
	diag.ub=inf;
	diag.delta_u=inf;
	diag.u_min=0;
	diag.u_max=0;
	diag.z_norm=0;
	diag.radius=node->anchor_radius;

	if(node->anchor_W<=0 || !isfinite(node->anchor_W))
		return diag;

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
		return diag;

	if(mu<u_min)
		mu=u_min;
	if(mu>u_max)
		mu=u_max;

	double anchor_L=expFromLogValue(logW-z2+mu);
	double anchor_U;

	if(u_max-u_min<epsilon)
	{
		anchor_U=expFromLogValue(logW-z2+u_min);
	}
	else
	{
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
	}

	if(anchor_U<anchor_L)
		anchor_U=anchor_L;

	diag.valid=isfinite(anchor_L) && anchor_L>=0.0 && anchor_L<=node->sumE && anchor_U>=anchor_L;
	diag.lb=anchor_L;
	diag.ub=anchor_U;
	diag.delta_u=(double)(u_max-u_min);
	diag.u_min=(double)u_min;
	diag.u_max=(double)u_max;
	diag.z_norm=(double)z_norm;

	return diag;
}

static double exactNodeContribution(kdAnchorAugNode*node,double*q,double**dataMatrix,int dim)
{
	double total=0;
	for(int i=0;i<(int)node->idList.size();i++)
	{
		int id=node->idList[i];
		double sq_dist=0;
		for(int d=0;d<dim;d++)
			sq_dist+=(q[d]-dataMatrix[id][d])*(q[d]-dataMatrix[id][d]);
		total+=exp(-sq_dist);
	}

	return total;
}

static double relGap(double lb,double ub)
{
	double denom=ub+lb;
	if(fabs(denom)<1e-300)
		return 0.0;

	return (ub-lb)/denom;
}

static void writeNodeDiagnostics(
	kdAnchorAugNode*node,
	double*q,
	double**dataMatrix,
	int dim,
	KDE_stat& stat,
	ofstream& out,
	int& next_id,
	int parent_id,
	int depth)
{
	int node_id=next_id++;
	double karl_lb=node->kdLinearAugNode::LB(q,dim,stat);
	double karl_ub=node->kdLinearAugNode::UB(q,dim,stat);
	double lb_dist=ell_MBR(q,node->boundary,dim);
	double ub_dist=u_MBR(q,node->boundary,dim);
	double delta_x=ub_dist*ub_dist-lb_dist*lb_dist;
	AnchorDiag anchor=computeAnchorDiag(node,q,dim);
	double exact=exactNodeContribution(node,q,dataMatrix,dim);

	double karl_gap=karl_ub-karl_lb;
	double anchor_gap=anchor.ub-anchor.lb;
	double karl_rel=relGap(karl_lb,karl_ub);
	double anchor_rel=relGap(anchor.lb,anchor.ub);
	double gap_ratio=(karl_gap>0 && anchor.valid) ? anchor_gap/karl_gap : inf;
	double rel_gap_ratio=(karl_rel>0 && anchor.valid) ? anchor_rel/karl_rel : inf;
	int karl_contains=(karl_lb<=exact+1e-8 && exact<=karl_ub+1e-8) ? 1 : 0;
	int anchor_contains=(anchor.valid && anchor.lb<=exact+1e-8 && exact<=anchor.ub+1e-8) ? 1 : 0;

	out<<node_id<<","
		<<parent_id<<","
		<<depth<<","
		<<node->idList.size()<<","
		<<(((int)node->idList.size()<=stat.n && node->childVector.empty()) ? 1 : 0)<<","
		<<node->sumE<<","
		<<delta_x<<","
		<<anchor.delta_u<<","
		<<anchor.z_norm<<","
		<<anchor.radius<<","
		<<karl_lb<<","
		<<karl_ub<<","
		<<karl_gap<<","
		<<karl_rel<<","
		<<anchor.valid<<","
		<<anchor.lb<<","
		<<anchor.ub<<","
		<<anchor_gap<<","
		<<anchor_rel<<","
		<<gap_ratio<<","
		<<rel_gap_ratio<<","
		<<exact<<","
		<<karl_contains<<","
		<<anchor_contains<<","
		<<((anchor.valid && anchor_gap<karl_gap) ? 1 : 0)<<","
		<<((anchor.valid && anchor_rel<karl_rel) ? 1 : 0)
		<<endl;

	for(int i=0;i<(int)node->childVector.size();i++)
		writeNodeDiagnostics((kdAnchorAugNode*)node->childVector[i],q,dataMatrix,dim,stat,out,next_id,node_id,depth+1);
}

int main(int argc,char**argv)
{
	if(argc<8)
	{
		cerr<<"Usage: "<<argv[0]<<" <query_file> <data_file> <output_csv> <leaf> <rel_error> <b> <query_index>"<<endl;
		return 1;
	}

	char*queryFile=argv[1];
	char*dataFile=argv[2];
	char*outputFile=argv[3];
	int leafCapacity=atoi(argv[4]);
	double rel_error=atof(argv[5]);
	double b=atof(argv[6]);
	int query_index=atoi(argv[7]);

	int dim;
	int qNum;
	double**queryMatrix;
	double**dataMatrix;
	KDE_stat stat={};
	stat.rel_error=rel_error;

	extract_FeatureVector(queryFile,qNum,dim,queryMatrix,stat);
	extract_FeatureVector(dataFile,stat.n,dim,dataMatrix,stat);

	if(query_index<0 || query_index>=qNum)
	{
		cerr<<"query_index out of range"<<endl;
		return 1;
	}

	preprocess_Data(dataMatrix,stat.n,dim,false,b);
	preprocess_Data(queryMatrix,qNum,dim,true,b);
	stat.qSquareNorm=sqNorm(queryMatrix[query_index],dim);

	kdTree tree(dim,dataMatrix,leafCapacity,stat);
	tree.rootNode=new kdAnchorAugNode();
	tree.build_kdTree(stat);
	tree.updateAugment((kdNode*)tree.rootNode);

	ofstream out(outputFile);
	if(!out.is_open())
	{
		cerr<<"Cannot open output file"<<endl;
		return 1;
	}

	out<<setprecision(17);
	out<<"node_id,parent_id,depth,count,is_leaf,sumE,delta_x,delta_u,z_norm,anchor_radius,"
		<<"karl_lb,karl_ub,karl_gap,karl_rel_gap,"
		<<"anchor_valid,anchor_lb,anchor_ub,anchor_gap,anchor_rel_gap,"
		<<"anchor_over_karl_gap,anchor_over_karl_rel_gap,exact,karl_contains,anchor_contains,"
		<<"anchor_abs_better,anchor_rel_better"<<endl;

	int next_id=0;
	writeNodeDiagnostics((kdAnchorAugNode*)tree.rootNode,queryMatrix[query_index],dataMatrix,dim,stat,out,next_id,-1,0);

	cout<<"Wrote "<<next_id<<" node diagnostics to "<<outputFile<<endl;
	return 0;
}
