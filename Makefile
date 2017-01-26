all:vtree_build vtree_demo
vtree_build: vtree_build.cpp
	g++  -O3 vtree_build.cpp  -lmetis -o vtree_build
vtree_demo: vtree_knn_demo.cpp
	g++ -O3 vtree_knn_demo.cpp -lmetis -o vtree_knn_demo



