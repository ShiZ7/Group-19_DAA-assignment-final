# Group-19_DAA-assignment-final



# Densest Subgraph via h-Clique Densities

This repository contains a C++ implementation of an exact algorithm to find the densest subgraph based on h-clique densities in an undirected graph.
It supports any h ≥ 2 and is tested on standard datasets like AS733 and Netscience.

## Description
Given an undirected graph G=(V,E), the goal is to find a subgraph S that maximizes the density defined by the number of h-cliques per node in S.
This project implements the exact algorithm based on a binary search framework combined with a flow-network construction to find the optimal densest subgraph.

The steps include:
* Finding all (h-1)-cliques.
* Constructing a specialized flow network.
* Running the Ford-Fulkerson algorithm to compute maximum flows and corresponding minimum cuts.
* Performing binary search over possible densities to find the best subgraph.


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







# Dense Subgraph Discovery via Clique Density


This project analyzes real-world network datasets (specifically AS733 and Netscience) to find dense subgraphs based on h-clique density. It uses clique enumeration, core decomposition, and Dinic's Max-Flow algorithm to process the data and identify densest subgraphs for different clique sizes (h values).


## Features
* Loads graph datasets in edge list format.
* Enumerates all (h-1)-cliques and h-cliques.
* Computes core decomposition for graph sparsification.
* Identifies the densest subgraphs based on h-clique density.
* Measures computation time and outputs density results for different values of h.
* Supports large sparse networks.

## Requirements
* C++17 or later
* A C++ compiler like g++ or clang++
* Standard C++ libraries (iostream, vector, queue, algorithm, etc.)

## Dataset 
* as733.txt — An AS-level Internet topology graph (available online).
* netscience.txt — Collaboration network among scientists (you'll need to similarly adapt the load_as733 function if using a different dataset format).

## Usage
1) Clone the repository
   git clone https://github.com/your-username/your-repo-name.git
   cd your-repo-name
2) Place the dataset (e.g., as733.txt) in the project directory.
3) Compile the code:
   g++ -std=c++17 -O2 -o dense_clique dense_clique.cpp
4) Run the program:
   ./dense_clique


The program will:
* Load the dataset.
* For h = 2, 3, ..., 6, find the densest subgraph based on h-clique density.
* Output the density results into a file called density_vs_h.txt.

5) Output:
   * Densities for each h are printed in the terminal.
   * The file density_vs_h.txt contains:

   h density
   2 0.0456
   3 0.0123
   4 0.0051
   ...
   where h is the clique size and density is the highest found density for subgraphs.

## File Structure
dense_clique.cpp      # Main source code
as733.txt             # Input graph file (example)
density_vs_h.txt      # Output file: density results

## Notes
* The code uses binary search to efficiently check for edge existence.
* For large graphs, clique enumeration can be computationally expensive for large h; you may want to adjust h accordingly.
* Currently, the code assumes that node indices are zero-based and contiguous (0, 1, 2, ...).

## Possible Improvements
* Extend to other datasets by modifying the graph loading function.
* Parallelize clique enumeration for large-scale graphs.
* Implement memory optimization for very large datasets.





### ID & Name: 

2022A7PS0047H
Shivansh Shanker Gupta

2020B4A71567H
T.Sai Sathwik

2022A7PS0041H
G Sri Vishwahitha

2022A7PS2014H
Snigdha Kaipa

2021B4A72488H
Sidharth Saxena


 ## Website Link

 https://shiz7.github.io/Group-19_DAA-assignment-final/




