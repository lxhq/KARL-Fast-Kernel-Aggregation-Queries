g++ -c init_KC.cpp -w -o init_KC.o -std=c++11
g++ -c SS_Binary.cpp -w -o SS_Binary.o -std=c++11
g++ -c Euclid_Bound.cpp -w -o Euclid_Bound.o -std=c++11
g++ -c m_tree.cpp -w -o m_tree.o -std=c++11
g++ -c kd_tree.cpp -w -o kd_tree.o -std=c++11
g++ -c Tree.cpp -w -o Tree.o -std=c++11
g++ -c GBF_KC_Binary.cpp -w -o GBF_KC_Binary.o -std=c++11
g++ -c svm.cpp -w -o svm.o -std=c++11
g++ -c my_svm.cpp -w -o my_svm.o -std=c++11

g++ main.cpp -O3 -o main init_KC.o SS_Binary.o Euclid_Bound.o m_tree.o kd_tree.o Tree.o GBF_KC_Binary.o svm.o my_svm.o -std=c++11

#####################################################################################################################################################################################################
#Type III-tau (Two-class SVM classification)
#	parameters in our main_KDC.cpp			(description of parameters)
#	int dim=atoi(argv[1]);				(dimensionality of the dataset)
#	int qNum=atoi(argv[2]);				(number of queries in query dataset)
#	char*querysetFileName=argv[3];			(query dataset name)
#	char*datasetFileName=argv[4];			(raw dataset name)
#	char*classResultName=argv[5];			(result file name)
#	int method=atoi(argv[6]);			(method=7 SCAN (with LibSVM), method=1 SOTA + kd-tree, method=11 KARL + kd-tree, method=5 SOTA + ball-tree, method=6 KARL + ball-tree)
#	int leafCapacity=atoi(argv[7]);			(leaf capacity of the tree-structure)
#	bulkLoad_TreeName=argv[8];			(Only used for ball-tree)

#Compare SCAN, SOTA and KARL using leaf capacity=80 in the online stage of 2-class SVM model
#SCAN
./main 22 10000 ./Datasets/ijcnn1/ijcnn1_sample.t ./Datasets/ijcnn1/ijcnn1_model_validate ./test_Result/ijcnn1_SCAN 7 80

#SOTA
./main 22 10000 ./Datasets/ijcnn1/ijcnn1_sample.t ./Datasets/ijcnn1/ijcnn1_model_validate ./test_Result/ijcnn1_SOTA_kd 1 80 #kd-tree
./main 22 10000 ./Datasets/ijcnn1/ijcnn1_sample.t ./Datasets/ijcnn1/ijcnn1_model_validate ./test_Result/ijcnn1_SOTA_ball 5 80 ./Datasets/ijcnn1/ijcnn1_SOTA_ball #ball-tree

#KARL
./main 22 10000 ./Datasets/ijcnn1/ijcnn1_sample.t ./Datasets/ijcnn1/ijcnn1_model_validate ./test_Result/ijcnn1_KARL_kd 11 80 #kd-tree
./main 22 10000 ./Datasets/ijcnn1/ijcnn1_sample.t ./Datasets/ijcnn1/ijcnn1_model_validate ./test_Result/ijcnn1_KARL_ball 6 80 ./Datasets/ijcnn1/ijcnn1_KARL_ball #ball-tree
#####################################################################################################################################################################################################
