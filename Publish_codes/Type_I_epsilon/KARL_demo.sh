g++ -c init_KAQ.cpp -w -o init_KAQ.o -std=c++11
g++ -c SS.cpp -w -o SS.o -std=c++11
g++ -c Euclid_Bound.cpp -w -o Euclid_Bound.o -std=c++11
g++ -c kd_tree.cpp -w -o kd_tree.o -std=c++11
g++ -c m_tree.cpp -w -o m_tree.o -std=c++11
g++ -c binaryTree.cpp -w -o binaryTree.o -std=c++11
g++ -c Validation.cpp -w -o Validation.o -std=c++11
g++ -c GBF_KAQ.cpp -w -o GBF_KAQ.o -std=c++11

g++ main.cpp -O3 -o main init_KAQ.o SS.o Euclid_Bound.o kd_tree.o m_tree.o binaryTree.o Validation.o GBF_KAQ.o

#####################################################################################################################################################################################################
#	parameters in our main.cpp			(description of parameters)
#	char*querysetFileName=argv[1];			(query dataset name)
#	char*datasetFileName=argv[2]; 			(raw dataset name)
#	char*resultFileName=argv[3]; 			(result file name)
#	int method=atoi(argv[4]); 			(method=0 SCAN, method=1 SOTA + kd-tree, method=3 KARL + kd-tree, method=11 SOTA + ball-tree, method=12 KARL + ball-tree)
#	int leafCapacity=atoi(argv[5]); 		(the leaf capacity) It can be 10,20,40,80,160,320 and 640
#	double rel_error=atof(argv[6]); 		(epsilon in our paper)
#	double b=atof(argv[7]); 			(the parameter b in KDE)

#Compare SCAN, SOTA and KARL using leaf capacity=80 with epsilon=0.2 in the online stage of kernel density estimation model

#SCAN
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_SCAN 0 80 0.2 10

#SOTA
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_SOTA_kd 1 80 0.2 10 #kd-tree
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_SOTA_ball 11 80 0.2 10 #ball-tree

#KARL
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd 3 80 0.2 10 #kd-tree
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball 12 80 0.2 10 #kd-tree


exit 
#####################################################################################################################################################################################################
#  	How to use auto-tuning?
#  	In HT_Sensor_qSet.dat, we only use 10000 queries, we sample 1000 queries in HT_Sensor dataset BUT NOT inside HT_Sensor_qSet.dat
#	We name the sample set is: HT_Sensor_sampleSet.dat (You can remove the above "exit" command once you finish this step)
#	Then, we run the following scripts in the offline stage:

./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_10_sample 3 10 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_20_sample 3 20 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_40_sample 3 40 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_80_sample 3 80 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_160_sample 3 160 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_320_sample 3 320 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_640_sample 3 640 0.2 10

./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_10_sample 12 10 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_20_sample 12 20 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_40_sample 12 40 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_80_sample 12 80 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_160_sample 12 160 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_320_sample 12 320 0.2 10
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_sampleSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_ball_640_sample 12 640 0.2 10

exit
#####################################################################################################################################################################################################
#After we obtain the best setting, we choose that parameters in the online stage. For example: if kd-tree and leaf capacity=160 can provide us the best throughput, we run the following script in the online stage (You can remove the above "exit" command once you obtain the best setting)
./main ./Datasets/HT_Sensor_UCI/HT_Sensor_qSet.dat ./Datasets/HT_Sensor_UCI/HT_Sensor_X.dat ./test_Result/HT_Sensor_UCI_method_KARL_kd_160 3 160 0.2 10

#For other models, we use the similar method for conducting the exeriments
