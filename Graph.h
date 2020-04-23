#ifndef GRAPH_HAS_BEEN_INCLUDED
#define GRAPH_HAS_BEEN_INCLUDED

#include "UTAdjacencyMatrix.h"

#include <queue>
#include <utility>
#include <vector>

namespace vspc
{

template <typename NodeType>
Graph
{
public:
    using Edge = std::pair<size_t, size_t>;

    Graph(const std::vector<NodeType>& nodes, const std::vector<Edge>& edges);

    void computeFundamentalCycles();

private:

    struct TreeNode {
        size_t    index;
        TreeNode* parent;
    };

    std::vector<NodeType> mNodes;
    UTAdjacencyMatrix mAdjMat;
};

// -----------------------------------------------------------------------------

template <typename NodeType>
Graph<NodeType>::Graph(const std::vector<NodeType&> nodes,
                   const std::vector<Edge>& edges)
    : mNodes(nodes)
    , mAdjMat(nodes.size())
{
    for (const auto& edge : edges) {
        mAdj.addEdge(edge.first, edge.second);
    }
}

template <typename NodeType>
Graph<NodeType>::computeFundamentalCycles()
{
    std::vector<TreeNode> spTree(mNodes.size());
    std::queue<size_t> nodeIdxQueue;

    // Arbitrarily start with the first node of the graph
    nodeIdxQueue.push(0);

    // Copy the adjacency matrix since we'll be removing edges
    UTAdjacencyMatrix adjMat = mAdjMat;

    // Initially, all tree nodes are its own parent
    for (size_t i = 0, n = mNodes.size(); i < n; ++i) {
        spTree[i].index  = i;
        spTree[i].parent = &(spTree[i]);
    }

    // TODO WIP
    // while (!nodeIdxQueue.empty()) {
    // }
}

} // namespace vspc

#endif // GRAPH_HAS_BEEN_INCLUDED
