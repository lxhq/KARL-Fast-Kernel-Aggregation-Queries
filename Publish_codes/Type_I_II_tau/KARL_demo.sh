g++ -c init_KC.cpp -w -o init_KC.o -std=c++11
g++ -c exp_LookUp.cpp -w -o exp_LookUp.o -std=c++11
g++ -c SS.cpp -w -o SS.o -std=c++11
g++ -c Euclid_Bound.cpp -w -o Euclid_Bound.o -std=c++11
g++ -c kd_tree.cpp -w -o kd_tree.o -std=c++11
g++ -c m_tree.cpp -w -o m_tree.o -std=c++11
g++ -c binaryTree.cpp -w -o binaryTree.o -std=c++11
g++ -c Oracle_Tree.cpp -w -o Oracle_Tree.o -std=c++11
g++ -c Tree.cpp -w -o Tree.o -std=c++11
g++ -c GBF_KC.cpp -w -o GBF_KC.o -std=c++11
g++ -c Sequential_Bound.cpp -w -o Sequential_Bound.o -std=c++11
g++ -c inverted_Index.cpp -w -o inverted_Index.o -std=c++11
g++ -c svm.cpp -w -o svm.o -std=c++11
g++ -c my_svm.cpp -w -o my_svm.o -std=c++11
g++ -c Incremental_Classification.cpp -w -o Incremental_Classification.o -std=c++11
g++ -c cache_order.cpp -w -o cache_order.o -std=c++11
g++ -c Reorder.cpp -w -o Reorder.o -std=c++11
g++ -c Stat_Estimation.cpp -w -o Stat_Estimation.o -std=c++11
g++ -c Oracle.cpp -w -o Oracle.o -std=c++11

g++ main_KDC.cpp -O3 -o main_KDC init_KC.o SS.o Euclid_Bound.o kd_tree.o m_tree.o Oracle_Tree.o Tree.o GBF_KC.o binaryTree.o Sequential_Bound.o inverted_Index.o exp_LookUp.o svm.o my_svm.o Incremental_Classification.o Reorder.o Stat_Estimation.o Oracle.o cache_order.o -std=c++11
g++ main.cpp -O3 -o main init_KC.o SS.o Euclid_Bound.o kd_tree.o m_tree.o Oracle_Tree.o Tree.o GBF_KC.o binaryTree.o Sequential_Bound.o inverted_Index.o exp_LookUp.o svm.o my_svm.o Incremental_Classification.o Reorder.o Stat_Estimation.o Oracle.o cache_order.o -std=c++11

#####################################################################################################################################################################################################
#Type I-tau (Kernel density classification)
#	parameters in our main_KDC.cpp			(description of parameters)
#	char*querysetFileName=argv[1];			(query dataset name)
#	char*datasetFileName=argv[2];			(raw dataset name)
#	char*classResultName=argv[3];			(result file name)
#	int method=atoi(argv[4]);			(method=0 SCAN, method=1 SOTA + kd-tree, method=4 KARL + kd-tree, method=35 SOTA + ball-tree, method=34 KARL + ball-tree)
#	int leafCapacity=atoi(argv[5]);			(leaf capacity of the tree-structure)
#	double rho=atof(argv[6]);			(threshold)
#	double b=atof(argv[7]);				(the parameter b in KDE)

#Compare SCAN, SOTA and KARL using leaf capacity=80 in the online stage of kernel density classification model

#SCAN
./main_KDC ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_SCAN_mu 0 80 106702 10

#SOTA
./main_KDC ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_SOTA_kd_mu 1 80 106702 10
./main_KDC ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_SOTA_ball_mu 35 80 106702 10

#KARL
./main_KDC ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_KARL_kd_mu 4 80 106702 10
./main_KDC ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_KARL_ball_mu 34 80 106702 10

#####################################################################################################################################################################################################
#Type II-tau (One-class SVM classification)
#	parameters in our main.cpp			(description of parameters)
#	int dim=atoi(argv[1]);				(dimensionality of the dataset)
#	int qNum=atoi(argv[2]);				(number of queries in query dataset)
#	char*querysetFileName=argv[3];			(query dataset name)
#	char*datasetFileName=argv[4];			(raw dataset name)
#	char*classResultName=argv[5];			(result file name)
#	int method=atoi(argv[6]);			(method=12 SCAN (with LibSVM), method=1 SOTA + kd-tree, method=4 KARL + kd-tree, method=35 SOTA + ball-tree, method=34 KARL + ball-tree)
#	int leafCapacity=atoi(argv[7]);			(leaf capacity of the tree-structure)

#Compare SCAN, SOTA and KARL using leaf capacity=80 in the online stage of 1-class SVM classification model

#SCAN
./main 41 10000 ./Datasets/NSL_KDD/NSLkdd_q_sample.dat ./Datasets/NSL_KDD/NSLkdd_model_26 ./test_Result/NSLkdd_method_SCAN 12 80 #SCAN with LibSVM

#SOTA
./main 41 10000 ./Datasets/NSL_KDD/NSLkdd_q_sample.dat ./Datasets/NSL_KDD/NSLkdd_model_26 ./test_Result/NSLkdd_method_SOTA_kd 1 80 #kd-tree
./main 41 10000 ./Datasets/NSL_KDD/NSLkdd_q_sample.dat ./Datasets/NSL_KDD/NSLkdd_model_26 ./test_Result/NSLkdd_method_SOTA_ball 35 80 #ball-tree

#KARL
./main 41 10000 ./Datasets/NSL_KDD/NSLkdd_q_sample.dat ./Datasets/NSL_KDD/NSLkdd_model_26 ./test_Result/NSLkdd_method_KARL_kd 4 80 #kd-tree
./main 41 10000 ./Datasets/NSL_KDD/NSLkdd_q_sample.dat ./Datasets/NSL_KDD/NSLkdd_model_26 ./test_Result/NSLkdd_method_KARL_ball 34 80 #ball-tree
#####################################################################################################################################################################################################
