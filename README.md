# V-Tree: Efficient kNN Search on Moving Objects with Road-Network Constraints

This project consists of implement of V-Tree, which designed to support kNN search on moving objects with road-network constraints. 

# [Dependence]
[METIS] is required before you run the program.it is used to partition the graph.
Thus, before compile our code, you must install METIS in your linux system.
METIS link & download: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
    ```
    get metics from
    http://glaros.dtc.umn.edu/gkhome/metis/metis/download
    ```
    
##[Example of Install METIS]
    ```
    wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
    gunzip metis-5.x.y.tar.gz
    tar -xvf metis-5.x.y.tar
    cd metics-5.x.y
    make config
    make
    make install
    ```

#[Compile]
--------------------------------------------
Simple build the file is make.

    ```
    make 
    ./vtree 
    ```
The manual compile is :
    ```
    g++ -O3 vtree.cpp -o vtree -lmetis
    ```

Beware the linking issue, in our case, we use "g++ ... -lmetis"
If it is not working, you can try "g++ ... -L/**/**/YOUR_METIS_LIB_PATH"
Make sure "metis.h" is in your default include(.h) directory.

#[Usage]
--------------------------------------------
1.Build a vtree index.
./vtree build_and_save ./input_edge_file ./output_file
2.Load VTree and run KNN test use 
./vtree load_and_run ./VTree_file $K(int) $car_percent(int) $change_percent(int)'
 k means the number of k neighborhoods 
 `car_percent` is the number of running vehicles on the road network
 `change_percent` means the percent of vehicles change to the other vertex between each query.

Some annotations were written among the code.


##[Input file format]
-------------------------------------------
The edge file is consist of two pars. First line is the overall information of the
 graph. The other line is the detail edge information.
    ```
    1089933 2545844    //The first line of the input file shows the number of vertices and edges.
    0 1 655            //The first row is the origin, the second row is destination, the third row is the weight of the edge.
    0 80 1393
    0 81 1426
    1 0 655
    1 2 1157
    1 6 5413
    2 1 1157
    2 3 2394
    2 7 3857
    ...
    ```
The VTree input is directed graph. If an undirected graph is used, it also support that. 
!The parameter RevE needs to be true when the undirected graph is used.
!If the number of graph is not 
In our code, we did not assert the input graph is connected graph
But connected graph must be guaranteed before METIS partition the graph
Hence, it is suggested you have to pre-process the input road network for your own dataset

##[Some useful parameter of the source code]
-----
Partition_Part       // the fanout of the vtree
Naive_Split_Limit    // the max number of the leaf node, which is τ+1 in the article.
NWFIX                // if the number of edge starts from 0 Set `NWFIX = true`, otherwise set it to false.

-----

For better understanding of our code, we provide example(New York Road Network dataset)
File use: (Note the file input format)
        nw.edge (graph edge file)


# Licence

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

