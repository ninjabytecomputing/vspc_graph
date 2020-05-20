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

    // Find SSSR and return the vector of cycles
    const auto& cycles = op.run();

    {
        std::ofstream fout(outputName);
        for (size_t i = 0, n = cycles.size(); i < n; ++i) {
            fout << convertGraphToPath(cycles[i]).strAsCSV() << "\n";
        }
    }

    return 0;
}

