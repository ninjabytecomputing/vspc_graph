#ifndef UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED
#define UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>


namespace vspc
{

class UndirectedGraph
{
public:
    using NodeType = uint32_t;
    using Edge     = std::pair<NodeType, NodeType>;

    using iterator       = std::map<NodeType, std::set<NodeType>>::iterator;
    using const_iterator = std::map<NodeType, std::set<NodeType>>::const_iterator;

    UndirectedGraph();
    explicit UndirectedGraph(const std::vector<Edge>& edges);

    void addNode(NodeType i);
    void removeNode(NodeType i);

    void addEdge(NodeType i, NodeType j);
    void removeEdge(NodeType i, NodeType j);

    /// @return @c True if nodes @a i and @j are connected.
    bool hasEdge(NodeType i, NodeType j) const;

    /// Remove all disconnected nodes.
    UndirectedGraph& prune();

    /// @return Set of indices of the nodes.
    std::set<NodeType> getNodes() const;
    /// @return Set of indices of the nodes that node @a i is connected to.
    std::set<NodeType> getConnections(NodeType i) const;

    /// XOR operation that overwrites this undirected graph.
    UndirectedGraph& operator^=(const UndirectedGraph& rhs);

    /// Swap undirected graphs.
    void swap(UndirectedGraph& other);

    iterator       begin()        { return mConnectivity.begin(); }
    const_iterator begin()  const { return mConnectivity.begin(); }
    const_iterator cbegin() const { return mConnectivity.cbegin(); }

    iterator       end()        { return mConnectivity.end(); }
    const_iterator end()  const { return mConnectivity.end(); }
    const_iterator cend() const { return mConnectivity.cend(); }

    /// @return Largest node index.
    NodeType maxNode() const { return mConnectivity.crbegin()->first; }
    /// @return Smallest node index.
    NodeType minNode() const { return mConnectivity.cbegin()->first; }
    /// @return Number of nodes in the graph.
    size_t numNodes() const { return mConnectivity.size(); }
    /// @return Number of edges in the graph.
    size_t numEdges() const { return mNumEdges; }
    /// @return Smallest node index.
    NodeType minNode() const { return mConnectivity.cbegin()->first; }
    /// @return Largest node index.
    NodeType maxNode() const { return mConnectivity.crbegin()->first; }

    /// @return Human-readable information about this undirected graph.
    std::string str() const;

private:
    size_t                                 mNumEdges;
    std::map<NodeType, std::set<NodeType>> mConnectivity;
};

// -----------------------------------------------------------------------------

UndirectedGraph operator^(UndirectedGraph lhs, const UndirectedGraph& rhs);

std::ostream& operator<<(std::ostream& os, const UndirectedGraph& obj);

// -----------------------------------------------------------------------------

UndirectedGraph::UndirectedGraph() : mNumEdges(0)
{
    // empty
}

UndirectedGraph::UndirectedGraph(const std::vector<Edge>& edges) : mNumEdges(0)
{
    for (const Edge& edge : edges) {
        addEdge(edge.first, edge.second);
    }
}

void
UndirectedGraph::addNode(NodeType i)
{
    mConnectivity[i];
}

void
UndirectedGraph::removeNode(NodeType i)
{
    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        mNumEdges -= it1->second.size();
        const auto itEnd = mConnectivity.erase(it1);
        for (auto it = mConnectivity.begin(); it != itEnd; ++it) {
            // Calling std::set::erase() on the key_type will return
            // the number of elements removed. Note that this value
            // should always be at most 1.
            mNumEdges -= it->second.erase(i);
        }
    }
}

void
UndirectedGraph::addEdge(NodeType i, NodeType j)
{
    // Don't insert self-loops
    if (i == j) return;

    // We only store the edge where i < j
    if (j < i) std::swap(i, j);

    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        const auto it2 = it1->second.find(j);
        if (it2 == it1->second.end()) {
            mConnectivity[i].insert(j);
            mConnectivity[j];
            ++mNumEdges;
        }
    } else {
        mConnectivity[i].insert(j);
        mConnectivity[j];
        ++mNumEdges;
    }
}

void
UndirectedGraph::removeEdge(NodeType i, NodeType j)
{
    // Don't insert self-loops
    if (i == j) return;

    // We only store the edge where i < j
    if (j < i) std::swap(i, j);

    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        const auto it2 = it1->second.find(j);
        if (it2 != it1->second.end()) {
            it1->second.erase(it2);
            --mNumEdges;
        }
    }
}

bool
UndirectedGraph::hasEdge(NodeType i, NodeType j) const
{
    // No self-loops
    if (i == j) return false;

    // We only store the edge where i < j
    if (j < i) std::swap(i, j);

    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        const auto it2 = it1->second.find(j);
        if (it2 != it1->second.end()) {
            return true;
        }
    }
    return false;
}

UndirectedGraph&
UndirectedGraph::prune()
{
    auto it = mConnectivity.begin();
    const auto itEnd = mConnectivity.end();
    while (it != itEnd) {
        // Aliases for convenience
        const NodeType&           node        = it->first;
        const std::set<NodeType>& connections = it->second;

        bool found = false;
        if (connections.empty()) {
            // Reverse lookup to see if this node is connected with any nodes
            // whose index are less than it.
            for (auto j = mConnectivity.begin(); j != it; ++j) {
                // Alias
                const std::set<NodeType>& s = j->second;
                if (s.find(node) != s.end()) {
                    found = true;
                    break;
                }
            }

            // If no connections were found, then remove the current node.
            if (!found) {
                it = mConnectivity.erase(it);
            } else {
                ++it;
            }
        } else {
            ++it;
        }
    }
    return *this;
}

std::set<UndirectedGraph::NodeType>
UndirectedGraph::getNodes() const
{
    std::set<NodeType> out;
    for (const auto& e : mConnectivity) {
        out.insert(e.first);
    }
    return out;
}

std::set<UndirectedGraph::NodeType>
UndirectedGraph::getConnections(NodeType i) const
{
    const auto it = mConnectivity.find(i);
    if (it != mConnectivity.end()) {
        return it->second;
    } else {
        return {};
    }
}

UndirectedGraph&
UndirectedGraph::operator^=(const UndirectedGraph& rhs)
{
    mNumEdges = 0;

    std::set<NodeType> rhsNodes = rhs.getNodes();

    for (auto& [myNode, myConnections] : mConnectivity) {
        const auto it1 = rhs.mConnectivity.find(myNode);
        if (it1 != rhs.mConnectivity.end()) {

            rhsNodes.erase(myNode);

            std::set<NodeType> symDiff;
            std::set_symmetric_difference(
                    myConnections.begin(), myConnections.end(),
                    it1->second.begin(), it1->second.end(),
                    std::inserter(symDiff, symDiff.begin()));

            mConnectivity[myNode] = symDiff;
            mNumEdges += symDiff.size();

        } else {
            // Do nothing, but count edges
            mNumEdges += mConnectivity[myNode].size();
        }
    }

    // RHS graph might still have nodes that this graph didn't have,
    // so just add in the edges.
    for (NodeType n : rhsNodes) {
        mConnectivity[n] = rhs.mConnectivity.at(n);
        mNumEdges += mConnectivity.at(n).size();
    }

    return *this;
}

void
UndirectedGraph::swap(UndirectedGraph& other)
{
    std::swap(mNumEdges, other.mNumEdges);
    std::swap(mConnectivity, other.mConnectivity);
}

std::string
UndirectedGraph::str() const
{
    // This method isn't trivial because this class only explicitly stores
    // edges via the smaller index of the two nodes that are connected.
    // Hence, in order to accurately display the connectivity of the entire
    // graph, some post-processing needs to happen to grab all relevant
    // information. For example, if nodes 1 and 3 are connected, then
    // it's trivial to list that node 1 is connected with node 3. However,
    // displaying the reverse information requires a reverse lookup.

    std::ostringstream os;
    os << "Number of nodes: " << numNodes() << "\n";
    os << "Number of edges: " << numEdges() << "\n";
    for (auto it = mConnectivity.begin(), itEnd = mConnectivity.end(); it != itEnd; ++it) {
        // Aliases for convenience
        const NodeType&           node        = it->first;
        const std::set<NodeType>& connections = it->second;

        // // Reverse lookup to see if this node is connected with any nodes
        // // whose index are less than it.
        // std::set<NodeType> tmp;
        // for (auto j = mConnectivity.begin(); j != it; ++j) {
        //     // Alias
        //     const std::set<NodeType>& s = j->second;
        //     if (s.find(node) != s.end()) {
        //         tmp.insert(j->first);
        //     }
        // }

        // // Merge what's found with the set of saved connections.
        // std::set<NodeType> allConnections;
        // std::set_union(tmp.begin(), tmp.end(),
        //                connections.begin(), connections.end(),
        //                std::inserter(allConnections, allConnections.begin()));

        // Write to stream
        os << "  Node: " << node;
        if (!connections.empty()) {
            os << " ---> ";
            for (NodeType n : connections) {
                os << n << " ";
            }
        }
        os << "\n";
    }
    return os.str();
}

UndirectedGraph
operator^(UndirectedGraph lhs, const UndirectedGraph& rhs)
{
    return lhs ^= rhs;
}

std::ostream&
operator<<(std::ostream& os, const UndirectedGraph& obj)
{
    return os << obj.str();
}

} // namespace vspc

#endif // UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED
