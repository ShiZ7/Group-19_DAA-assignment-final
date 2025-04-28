# Group-19_DAA-assignment-final



## 1. Requirements
* A C++17 compatible compiler (e.g., g++, clang++)
* Linux, MacOS, or Windows with a terminal
* (Optional) make for easier building

## 2. Compilation2. Compilation
g++ -std=c++17 -O2 -o densest_subgraph main.cpp

## 3. Running
./densest_subgraph <graph_file> <h>  , where
* <graph_file>: Path to the input graph file.
* <h>: The size of cliques to consider (h ≥ 2).


## Input Format
n m
u1 v1
u2 v2
...
um vm

Where:
* n = number of nodes
* m = number of edges
* Each line after gives an undirected edge between nodes u and v.
* Nodes are indexed from 0 to n-1.


## Output

After execution, the program prints:
* Number of nodes and edges in the loaded graph.
* Details about the search process (binary search steps).
* The densest subgraph found:
* Number of nodes
* Number of h-cliques
* Clique-density
* List of nodes in the subgraph
* Total execution time in milliseconds.

## Datasets Used
* AS733
* Netscience

## Performance Notes
* The exact algorithm guarantees finding the optimal subgraph, but can be computationally heavy for large graphs or large values of h.
* Works efficiently for moderate-sized graphs (up to a few thousand nodes) and small h values (2, 3, or 4)










