#include "Graph.h"

#include <iostream>

int main() {
    vspc::Graph g;
    g.addEdge(1, 2);
    g.addEdge(2, 1);    // Should not be included
    g.addEdge(1, 2);    // Should not be included
    g.addEdge(1, 3);
    g.addEdge(2, 3);
    g.addEdge(1, 4);
    g.addEdge(1, 10);

    std::cout << g << "\n";

    std::cout << g.isConnected(1, 2) << std::endl;
    std::cout << g.isConnected(2, 1) << std::endl;
    std::cout << g.isConnected(1, 11) << std::endl;
}
