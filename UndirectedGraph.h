#ifndef UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED
#define UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace vspc
{

class UndirectedGraph
{
public:
    using Edge = std::pair<int, int>;

    UndirectedGraph();
    explicit UndirectedGraph(const std::vector<Edge>& edges);

    void addNode(int i);
    void removeNode(int i);

    void addEdge(int i, int j);
    void removeEdge(int i, int j);

    bool isConnected(int i, int j) const;

    UndirectedGraph& operator^=(const UndirectedGraph& rhs);

    size_t numNodes() const { return mConnectivity.size(); }
    size_t numEdges() const { return mNumEdges; }

    std::string str() const;

private:
    size_t                       mNumEdges;
    std::map<int, std::set<int>> mConnectivity;
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
UndirectedGraph::addNode(int i)
{
    mConnectivity[i];
}

void
UndirectedGraph::removeNode(int i)
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
UndirectedGraph::addEdge(int i, int j)
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
UndirectedGraph::removeEdge(int i, int j)
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
UndirectedGraph::isConnected(int i, int j) const
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
UndirectedGraph::operator^=(const UndirectedGraph& rhs)
{
    mNumEdges = 0;
    for (auto&& [myNode, myConnections] : mConnectivity) {
        const auto it1 = rhs.mConnectivity.find(myNode);
        if (it1 != rhs.mConnectivity.end()) {
            std::set<int> symDiff;
            std::set_symmetric_difference(
                    myConnections.begin(), myConnections.end(),
                    it1->second.begin(), it1->second.end(),
                    std::inserter(symDiff, symDiff.begin()));

            mConnectivity[myNode] = symDiff;
            mNumEdges += symDiff.size();
        }
    }
    return *this;
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
    os << "Undirected Graph:\n";
    os << "  Number of nodes: " << numNodes() << "\n";
    os << "  Number of edges: " << numEdges() << "\n";
    for (auto it = mConnectivity.begin(), itEnd = mConnectivity.end(); it != itEnd; ++it) {
        // Aliases for convenience
        const int&           node        = it->first;
        const std::set<int>& connections = it->second;

        os << "  Node: " << node;
        if (!connections.empty()) {
            os << " ---> ";
            for (auto&& n : connections) {
                os << n << " ";
            }
        } else {
            // Reverse lookup to see if this node is connected with any nodes
            // whose index are less than it.
            std::set<int> tmp;
            for (auto j = mConnectivity.begin(); j != it; ++j) {
                // Alias
                const auto& s = j->second;
                if (s.find(node) != s.end()) {
                    tmp.insert(j->first);
                }
            }
            if (!tmp.empty()) {
                os << " ---> ";
                for (auto&& n : tmp) {
                    os << n << " ";
                }
            }
        }
        os << "\n";
    }
    return os.str();
}

// -----------------------------------------------------------------------------

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
