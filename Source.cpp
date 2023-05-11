#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include <chrono>
#include <fstream>
#include <limits>
using namespace std::chrono;
std::fstream out("results.txt", std::fstream::in | std::fstream::out | std::fstream::app);
const int INF = std::numeric_limits<int>::max();

struct Edge {
    int src, dest, weight;
};

int findMinKeyVertex(const std::vector<int>& key, const std::vector<bool>& mstSet, int numNodes) {
    int minKey = INF;
    int minIndex = -1;

    for (int v = 0; v < numNodes; ++v) {
        if (!mstSet[v] && key[v] < minKey) {
            minKey = key[v];
            minIndex = v;
        }
    }

    return minIndex;
}

// Function to perform Prim's algorithm on a graph using adjacency matrix representation
void primMST(const std::vector<std::vector<int>>& graph, int numNodes) {
    std::vector<int> key(numNodes);     // Array to store the key value of each node
    std::vector<bool> mstSet(numNodes); // Array to track if a node is already included in the MST

    // Initialize all keys as infinity and MST set as false
    for (int v = 0; v < numNodes; ++v) {
        key[v] = INF;
        mstSet[v] = false;
    }

    // The first vertex is always the root of the MST
    key[0] = 0;

    // Construct the MST with (numNodes - 1) edges
    for (int count = 0; count < numNodes - 1; ++count) {
        // Find the vertex with the minimum key value from the set of vertices not yet included in MST
        int u = findMinKeyVertex(key, mstSet, numNodes);

        // Mark the selected vertex as visited
        mstSet[u] = true;

        // Update the key values of the adjacent vertices of the selected vertex
        for (int v = 0; v < numNodes; ++v) {
            if (graph[u][v] != 0 && !mstSet[v] && graph[u][v] < key[v]) {
                key[v] = graph[u][v];
            }
        }
    }

    // Calculate and return the cost of the minimum spanning tree
    int cost = 0;
    for (int v = 0; v < numNodes; ++v) {
        cost += key[v];
    }

    std::cout << "Final Cost of Minimum Spanning Tree(PRIM): " << cost << std::endl;
}

// Function to compare edges based on their weights (for sorting)
bool compareEdges(const Edge& e1, const Edge& e2) {
    return e1.weight < e2.weight;
}

// Function to find the parent of a node in the disjoint set
int findParent(std::vector<int>& parent, int node) {
    if (parent[node] == node)
        return node;
    return findParent(parent, parent[node]);
}

// Function to perform union operation of two sets (by rank)
void unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
    int xParent = findParent(parent, x);
    int yParent = findParent(parent, y);

    if (rank[xParent] < rank[yParent])
        parent[xParent] = yParent;
    else if (rank[xParent] > rank[yParent])
        parent[yParent] = xParent;
    else {
        parent[yParent] = xParent;
        rank[xParent]++;
    }
}

// Function to perform Kruskal's algorithm on a graph using adjacency matrix representation
void kruskalMST(const std::vector<std::vector<int>>& graph, int numNodes) {
    // Create a vector to store the edges
    std::vector<Edge> edges;

    // Traverse the upper triangular portion of the adjacency matrix and store the edges
    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            if (graph[i][j] != 0) {
                edges.push_back({ i, j, graph[i][j] });
            }
        }
    }

    // Sort the edges in ascending order based on their weights
    std::sort(edges.begin(), edges.end(), compareEdges);

    // Create a vector to store the parent of each node in the disjoint set
    std::vector<int> parent(numNodes);

    // Initialize the parent vector such that each node is its own parent initially
    for (int i = 0; i < numNodes; ++i) {
        parent[i] = i;
    }

    // Create a vector to store the rank of each node in the disjoint set (for union by rank)
    std::vector<int> rank(numNodes, 0);

    // Create a vector to store the edges of the minimum spanning tree
    std::vector<Edge> mst;

    // Traverse through all the edges in the sorted order
    for (const auto& edge : edges) {
        int srcParent = findParent(parent, edge.src);
        int destParent = findParent(parent, edge.dest);

        // If including this edge doesn't form a cycle, add it to the minimum spanning tree
        if (srcParent != destParent) {
            mst.push_back(edge);
            unionSets(parent, rank, srcParent, destParent);
        }
    }

    // Calculate the final cost of the minimum spanning tree
    int cost = 0;
    for (const auto& edge : mst) {
        cost += edge.weight;
    }

    // Display the final cost of the minimum spanning tree
    std::cout << "Final Cost of Minimum Spanning Tree(KRUSKAL): " << cost << std::endl;
}

// Function to generate a random weighted graph with a specified number of nodes
std::vector<std::vector<int>> generateRandomWeightedGraph(int numNodes) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 10);

    std::vector<std::vector<int>> graph(numNodes, std::vector<int>(numNodes, INF));

    // Generate random weights for the graph
    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            int weight = dist(gen);
            graph[i][j] = weight;
            graph[j][i] = weight;
        }
    }

    // Set diagonal elements to 0
    for (int i = 0; i < numNodes; ++i) {
        graph[i][i] = 0;
    }

    return graph;
}

// Function to print the adjacency matrix of a graph
void printGraph(const std::vector<std::vector<int>>& graph) {
    for (const auto& row : graph) {
        for (int val : row) {
            if (val == INF)
                std::cout << "INF ";
            else
                std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void outputMetrics(char f, std::vector<std::vector<int>>graph, int vert)
{
    std::string s;
    auto search_start = high_resolution_clock::now();
    switch (f)
    {
    case '1':
        primMST(graph, vert);
        s = "Prim: ";
        break;
    case '2':
        kruskalMST(graph, vert);
        s = "Kruskal: ";
        break;
    }
    auto search_stop = high_resolution_clock::now();
    auto duration_search = duration_cast<milliseconds>(search_stop - search_start);
    out  << duration_search.count() / 1e3f << std::endl;
}


int main() {
    int numNodes;
    
    //std::cout << "Random Weighted Graph (Adjacency Matrix):" << std::endl;
    //printGraph(graph);

    //kruskalMST(graph, numNodes);
    //primMST(graph1, numNodes);
    for (int i = 50; i <= 5500; i += 50)
        //{
            //numNodes = 5500;
            //std::vector<std::vector<int>> graph = generateRandomWeightedGraph(numNodes);
            //std::vector<std::vector<int>> graph1 = graph;

            //outputMetrics('1', graph, numNodes);
            //outputMetrics('2', graph, numNodes);
        std::cout << i << std::endl;

    //}
    return 0;
}
