

In this package, we used the source code from the paper of 
Ken C. K. Lee, Wang-Chien Lee, and Baihua Zheng. Fast object search on road networks. In EDBT, pages 1018–1029, 2009.
The authors of the paper have all the right of their source code. We just do some modification to evaluate the performance of kNN searching. 
We added some files to evaluate the performance.
hiernn_vtree.cc
hiernn_vtree_density.cc
hiernn_vtree_dif_dist_test.cc

To compile the file with 
$make
The bin file store the program used to evaluate the performance.

For example:
1 step . Generate the index for x road network
./hiergraphloader -n x.cnode -e x.cedge  -t 4 -l 8 -h xx_roadtree_t4_l8.idx  -v 
2 step . Load the index and test different performance.

./hiernn_vtree -h xx_roadtree_t4_l8.idx -x testFile -k ${k}

More detail can be found in the source code.
