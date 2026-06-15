#include "weak_auto.h"

int weak_auto(double**queryMatrix,double**dataMatrix,double*alphaArray,svm_node**svm_qMatrix,svm_model*& model,int sNum,int dim,double sum_alpha,double*a_G,double S_G,double*center,double radius,SVM_stat& stat,int leafCapacity,kdTree_adv& kd_Tree_adv)
{
	double seq_Time;
	double kd_Time;

	#ifdef C_PLUSPLUS11_CLOCK
		auto start_seq_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_seq_s=clock();
	#endif

	for(int s=0;s<sNum;s++)
		sequential_bound(queryMatrix[s],dataMatrix,alphaArray,svm_qMatrix[s],model,dim,sum_alpha,a_G,S_G,center,radius,stat,9);

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_seq_s=chrono::high_resolution_clock::now();
		seq_Time=(chrono::duration_cast<chrono::nanoseconds>(end_seq_s-start_seq_s).count())/1000000000.0;
	#else
		clock_t end_seq_s=clock();
		seq_Time=((double)(end_seq_s-start_seq_s))/CLOCKS_PER_SEC;
	#endif


	#ifdef C_PLUSPLUS11_CLOCK
		auto start_kd_s=chrono::high_resolution_clock::now();
	#else
		clock_t start_kd_s=clock();
	#endif

	for(int s=sNum;s<2*sNum;s++)
	{
		stat.qSquareNorm=computeSqNorm(queryMatrix[s],dim);
		GBF_iter(queryMatrix[s],svm_qMatrix[s],kd_Tree_adv,dim,stat);
	}

	#ifdef C_PLUSPLUS11_CLOCK
		auto end_kd_s=chrono::high_resolution_clock::now();
		kd_Time=(chrono::duration_cast<chrono::nanoseconds>(end_seq_s-start_seq_s).count())/1000000000.0;
	#else
		clock_t end_kd_s=clock();
		kd_Time=((double)(end_seq_s-start_seq_s))/CLOCKS_PER_SEC;
	#endif

	if(seq_Time<kd_Time)
		return 9;
	else
		return 36;
}