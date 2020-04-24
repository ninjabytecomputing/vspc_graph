#include "UndirectedGraph.h"

#include <iostream>

// #define FIRST_TEST
#define SECOND_TEST

int main() {

#ifdef FIRST_TEST
    vspc::UndirectedGraph g;
    g.addEdge(1, 2);
    g.addEdge(2, 1);    // Should not be included
    g.addEdge(1, 2);    // Should not be included
    g.addEdge(3, 1);
    g.addEdge(2, 3);
    g.addEdge(1, 4);
    g.addEdge(10, 1);
    g.addEdge(1, 1);    // Should not be included
    g.addNode(12);
    g.addNode(13);
    g.addEdge(13, 2);
    g.addNode(1);

    std::cout << g << std::endl;

    g.removeEdge(2, 1);
    g.removeEdge(1, 2);    // Should not do anything
    g.removeEdge(1, 3);
    g.removeEdge(2, 3);
    g.removeEdge(1, 4);
    g.removeEdge(12, 1);
    // g.removeNode(13);

    std::cout << g << std::endl;
#endif

#ifdef SECOND_TEST
    vspc::UndirectedGraph g1, g2;
    g1.addEdge(1, 2);
    g1.addEdge(3, 1);

    g2.addEdge(3, 2);
    g2.addEdge(2, 1);

    vspc::UndirectedGraph out = g1 ^ g2;

    std::cout << out << std::endl;
#endif
}
