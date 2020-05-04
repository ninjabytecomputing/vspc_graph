#include <SSSR.h>
#include <UndirectedGraph.h>
#include <UpperTriangularMatrix.h>

#include <iostream>

// #define TWO_RINGS
// #define PAPER_EXAMPLE
#define THREE_RINGS

int main() {
    vspc::UndirectedGraph graph;
    // Two triangles
#ifdef TWO_RINGS
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 0);
    graph.addEdge(1, 3);
#endif

#ifdef PAPER_EXAMPLE
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(4, 5);
    graph.addEdge(5, 6);
    graph.addEdge(6, 0);
    graph.addEdge(1, 5);
#endif

#ifdef THREE_RINGS
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(0, 4);

    graph.addEdge(1000, 1);
    graph.addEdge(2, 50);
    graph.addEdge(1000, 50);
    graph.addEdge(1, 2);

    graph.addEdge(3, 50);
    graph.addEdge(2, 3);
    graph.addEdge(2, 50);
#endif

    std::cout << "#############################################" << std::endl;
    std::cout << "# Input graph" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << graph << std::endl;

    vspc::SSSR op(graph);
    const auto& sssrSet = op.run();

    std::cout << "#############################################" << std::endl;
    std::cout << "# SSSR" << std::endl;
    std::cout << "#############################################" << std::endl;
    for (const vspc::UndirectedGraph& cycle : sssrSet) {
        std::cout << cycle << std::endl;
    }
}
