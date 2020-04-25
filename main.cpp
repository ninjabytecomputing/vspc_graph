#include "UndirectedGraph.h"

#include <iostream>

// #define FIRST_TEST
#define SECOND_TEST
// #define THIRD_TEST

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
    g1.addEdge(0, 1);
    g1.addEdge(1, 2);
    g1.addEdge(2, 6);

    g2.addEdge(0, 1);
    g2.addEdge(1, 5);
    g2.addEdge(5, 6);

    vspc::UndirectedGraph out = g1 ^ g2;
    out.prune();

    std::cout << g1 << std::endl;
    std::cout << g2 << std::endl;
    std::cout << out << std::endl;
#endif

#ifdef THIRD_TEST
    vspc::UndirectedGraph graph;
    // Case 1
    // graph.addEdge(0, 1);
    // graph.addEdge(0, 2);
    // graph.addEdge(0, 3);
    // graph.addEdge(1, 2);
    // graph.addEdge(2, 3);

    // Case 2
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(0, 4);

    graph.addEdge(5, 1);
    graph.addEdge(2, 6);
    graph.addEdge(5, 6);
    graph.addEdge(1, 2);

    // graph.addEdge(3, 6);
    // graph.addEdge(2, 3);
    // graph.addEdge(2, 6);

    std::cout << "Original graph" << std::endl;
    std::cout << graph << std::endl;

    vspc::FundamentalCyclesOp op(graph);
    op.compute();
    const std::vector<vspc::UndirectedGraph>& cycles = op.getFundamentalCycles();

    for (const auto& c : cycles) {
        std::cout << "Fundamental cycle:" << std::endl;
        std::cout << c << std::endl;
    }
#endif

}
