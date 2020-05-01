#include <SSSR.h>
#include <UndirectedGraph.h>
#include <UpperTriangularMatrix.h>

#include <iostream>

#define TWO_RINGS
// #define PAPER_EXAMPLE
// #define THREE_RINGS

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

    graph.addEdge(5, 1);
    graph.addEdge(2, 6);
    graph.addEdge(5, 6);
    graph.addEdge(1, 2);

    graph.addEdge(3, 6);
    graph.addEdge(2, 3);
    graph.addEdge(2, 6);
#endif

    std::cout << graph << std::endl;

    std::cout << "Initializing SSSR" << std::endl;
    vspc::SSSR sssrOp(graph);
    std::cout << "Initializing PID matrices" << std::endl;
    sssrOp.initializePID();
    std::cout << "Detecting candidates" << std::endl;
    sssrOp.makeCandidateSet();
}
