#include "Library.h"

int myrandom (int i) { return rand()%i;}

/*int main(int argc,char**argv)
{
	int num_of_generation=atoi(argv[1]);
	int total_sv=atoi(argv[2]);
	char*sample_idFileName=argv[3];

	//int num_of_generation=200;
	//int total_sv=97290;
	//char*sample_idFileName=(char*)"KDDCup_sample_id.txt";
	
	vector<int> sv_idVector;
	fstream sample_idFile;

	for(int i=0;i<total_sv;i++)
		sv_idVector.push_back(i);

	srand(time(0));
	random_shuffle(sv_idVector.begin(),sv_idVector.end(),myrandom);

	sample_idFile.open(sample_idFileName,ios::in|ios::out|ios::trunc);
	sample_idFile<<num_of_generation<<endl;
	
	for(int i=0;i<num_of_generation;i++)
		sample_idFile<<sv_idVector[i]<<endl;

	sample_idFile.close();
}*/