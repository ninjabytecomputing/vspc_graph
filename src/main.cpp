#include "CSVReader.h"
#include "UndirectedGraph.h"
#include "SSSR.h"

#include <fstream>
#include <string>

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

    std::cout << "+------------------------------------------------------+\n";
    std::cout << "|                      Input data                      |\n";
    std::cout << "+------------------------------------------------------+\n";
    std::cout << "  File: " << argv[1] << "\n";
    std::cout << "    Number of nodes : " << graph.numNodes() << "\n";
    std::cout << "    Node index range: "
              << "[" << graph.minNode() << ", " << graph.maxNode() << "]\n";
    std::cout << "    Number of edge  : " << graph.numEdges() << "\n\n";

    const auto& cycles = op.run();

    std::cout << "+------------------------------------------------------+\n";
    std::cout << "|                        Results                       |\n";
    std::cout << "+------------------------------------------------------+\n";
    std::cout << "  Found " << cycles.size() << " cycles\n";

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

