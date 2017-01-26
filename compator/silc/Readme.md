##SILC
 SILC can not support large road network, we only used it to evaluate to synetatic data. 
##To compile the silc 
make
The progress is also to generate the index first and then evaluate the performance.

build:
	g++ -O3 silc_build.cpp -o silc_build
query:
	g++ -O3 silc_knn.cpp -o silc_knn

##Usage example:
g++ -O3 silc_build.cpp –o silc_build

./silc_build cal.node cal.edge 1 cal.morton 1 

./silc_knn