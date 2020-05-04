#include "CSVReader.h"
#include "UndirectedGraph.h"
#include "SSSR.h"

#include <fstream>
#include <string>

vspc::SSSR::Path convertGraphToPath(vspc::UndirectedGraph graph);

int main(int argc, char* argv[]) {

    if (argc == 1) {
        std::cerr << "Need to specify a csv file containing graph data\n";
        return 1;
    }

    std::vector<vspc::UndirectedGraph::Edge> edgeData;

    {
        std::ifstream fin(argv[1]);
        if (fin.fail()) {
            std::cerr << "Invalid filename\n";
        }

        vspc::CSVIterator iter(fin);

        // Iterate past the metadata
        while (iter->isMetadata()) { ++iter; }

        for (; iter; ++iter) {
            const vspc::CSVRow& edge = *iter;
            edgeData.emplace_back(std::stoi(edge[0]), std::stoi(edge[1]));
        }
    }

    // Construct graph object
    vspc::UndirectedGraph graph(edgeData);
    // Construct SSSR object
    vspc::SSSR op(graph);

    std::cout << "+----------------------------------------------+\n";
    std::cout << "|                  Input data                  |\n";
    std::cout << "+----------------------------------------------+\n";
    std::cout << "  File: " << argv[1] << "\n";
    std::cout << "    Number of nodes : " << graph.numNodes() << "\n";
    std::cout << "    Node index range: "
              << "[" << graph.minNode() << ", " << graph.maxNode() << "]\n";
    std::cout << "    Number of edge  : " << graph.numEdges() << "\n";
    std::cout << "    Number of cycles: " << op.numTheoreticalCycles() << "\n\n";

    const auto& cycles = op.run();

    std::cout << "+----------------------------------------------+\n";
    std::cout << "|                    Results                   |\n";
    std::cout << "+----------------------------------------------+\n";
    std::cout << "  Found " << cycles.size() << "/"
                            << op.numTheoreticalCycles() << " cycle(s)\n";

    for (size_t i = 0, n = cycles.size(); i < n; ++i) {
        // Sanity check!
        if (cycles[i].numEdges() != cycles[i].numNodes()) {
            std::cout << "Cycle " << i << " failed the sanity check\n";
            std::cout << cycles[i] << std::endl;
        } else {
            // std::cout << cycles[i] << std::endl;
            std::cout << convertGraphToPath(cycles[i]) << "\n";
        }
    }
}

vspc::SSSR::Path convertGraphToPath(vspc::UndirectedGraph graph)
{
    using NodeType = vspc::UndirectedGraph::NodeType;
    using Path     = vspc::SSSR::Path;

    NodeType node = graph.minNode();
    std::set<NodeType> conn = graph.getConnections(node);

    Path path(node, *conn.cbegin());
    node = *conn.cbegin();

    while (graph.prune().numNodes() != 0) {
        conn = graph.getConnections(node);
        if (conn.empty()) {
            for (auto&& [n, c] : graph) {
                const auto it = c.find(node);
                if (it != c.end()) {
                    path.append(Path({node, n}));
                    graph.removeNode(node);
                    node = n;
                    break;
                }
            }
        } else {
            path.append(Path{node, *conn.cbegin()});
            graph.removeNode(node);
            node = *conn.cbegin();
        }
    }

    return path;
}
