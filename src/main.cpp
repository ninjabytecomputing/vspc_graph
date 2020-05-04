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

    for (const auto& e : edgeData) {
        std::cout << e.first << " " << e.second << std::endl;
    }
}
