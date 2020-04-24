#ifndef GRAPH_HAS_BEEN_INCLUDED
#define GRAPH_HAS_BEEN_INCLUDED

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace vspc
{

class Graph
{
public:
    using Edge = std::pair<int, int>;

    Graph();
    explicit Graph(const std::vector<Edge>& edges);

    // Default copy, move, and destructor are fine

    void addEdge(const int i, const int j);

    bool isConnected(const int i, const int j) const;

    size_t numNodes() const { return mConnectivity.size(); }
    size_t numEdges() const { return mNumEdges; }

    std::string str() const;

private:
    size_t                       mNumEdges;
    std::map<int, std::set<int>> mConnectivity;
};

std::ostream& operator<<(std::ostream& os, const Graph& obj);

// -----------------------------------------------------------------------------

Graph::Graph() : mNumEdges(0)
{
    // empty
}

Graph::Graph(const std::vector<Edge>& edges) : mNumEdges(0)
{
    for (const Edge& edge : edges) {
        addEdge(edge.first, edge.second);
    }
}

void
Graph::addEdge(const int i, const int j)
{
    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        const auto it2 = it1->second.find(j);
        if (it2 == it1->second.end()) {
            mConnectivity[i].insert(j);
            mConnectivity[j].insert(i);
            ++mNumEdges;
        }
    } else {
        mConnectivity[i].insert(j);
        mConnectivity[j].insert(i);
        ++mNumEdges;
    }
}

bool
Graph::isConnected(const int i, const int j) const
{
    const auto it1 = mConnectivity.find(i);
    if (it1 != mConnectivity.end()) {
        const auto it2 = it1->second.find(j);
        if (it2 != it1->second.end()) {
            return true;
        }
    }
    return false;
}

std::string
Graph::str() const
{
    std::ostringstream os;
    os << "Graph:\n";
    os << "  Number of nodes: " << numNodes() << "\n";
    os << "  Number of edges: " << numEdges() << "\n";
    for (auto&& [node, connections] : mConnectivity) {
        os << "  Node: " << node << " ---> ";
        for (auto&& n : connections) {
            os << n << " ";
        }
        os << "\n";
    }
    return os.str();
}

std::ostream&
operator<<(std::ostream& os, const Graph& obj)
{
    return os << obj.str();
}

} // namespace vspc

#endif // GRAPH_HAS_BEEN_INCLUDED
