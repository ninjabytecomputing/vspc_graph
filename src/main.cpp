#include "CSVReader.h"
#include "UndirectedGraph.h"
#include "SSSR.h"

#include <fstream>
#include <string>

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << "Need both input and output csv files\n";
        return 1;
    }

    std::string inputName(argv[1]);
    std::string outputName(argv[2]);
    std::vector<vspc::UndirectedGraph::Edge> edgeData;

    {
        std::ifstream fin(inputName);
        if (fin.fail()) {
            std::cerr << "Invalid filename: " << inputName << "\n";
            return 1;
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

    std::cout << "+------------------------------------------------------------------+\n";
    std::cout << "|                            Input data                            |\n";
    std::cout << "+------------------------------------------------------------------+\n";
    std::cout << "  Input file       : " << inputName << "\n";
    std::cout << "  Number of nodes  : " << graph.numNodes() << "\n";
    std::cout << "  Node index range : "
              << "[" << graph.minNode() << ", " << graph.maxNode() << "]\n";
    std::cout << "  Number of edges  : " << graph.numEdges() << "\n\n";

    // Find SSSR and return the vector of cycles
    const auto& cycles = op.run();

    std::cout << "+------------------------------------------------------------------+\n";
    std::cout << "|                              Results                             |\n";
    std::cout << "+------------------------------------------------------------------+\n";
    std::cout << "  Ouptut file      : " << outputName << "\n";
    std::cout << "  Number of cycles : " << cycles.size() << "\n";

    {
        std::ofstream fout(outputName);

        for (size_t i = 0, n = cycles.size(); i < n; ++i) {
            // Sanity check!
            if (cycles[i].numEdges() != cycles[i].numNodes()) {
                std::cout << "Cycle " << i << " failed the sanity check - skipping\n";
            } else {
                fout << convertGraphToPath(cycles[i]).strAsCSV() << "\n";
            }
        }
    }
}

