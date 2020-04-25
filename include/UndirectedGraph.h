#ifndef UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED
#define UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include <iostream>

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

    bool hasEdge(int i, int j) const;

    void prune();

    std::set<int> getNodes() const;
    std::set<int> getConnections(int i) const;

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
UndirectedGraph::hasEdge(int i, int j) const
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

void
UndirectedGraph::prune()
{
    auto it = mConnectivity.begin();
    const auto itEnd = mConnectivity.end();
    while (it != itEnd) {
        // Aliases for convenience
        const int&           node        = it->first;
        const std::set<int>& connections = it->second;

        if (connections.empty()) {
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

            if (tmp.empty()) {
                // Erase automatically moves the iterator to the
                // next element
                it = mConnectivity.erase(it);
            } else {
                ++it;
            }
        }
    }
}

std::set<int>
UndirectedGraph::getNodes() const
{
    std::set<int> out;
    for (const auto& e : mConnectivity) {
        out.insert(e.first);
    }
    return out;
}

std::set<int>
UndirectedGraph::getConnections(int i) const
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

    std::set<int> rhsNodes = rhs.getNodes();

    for (auto&& [myNode, myConnections] : mConnectivity) {
        const auto it1 = rhs.mConnectivity.find(myNode);
        if (it1 != rhs.mConnectivity.end()) {

            rhsNodes.erase(myNode);

            std::set<int> symDiff;
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
    for (int n : rhsNodes) {
        mConnectivity[n] = rhs.mConnectivity.at(n);
        mNumEdges += mConnectivity.at(n).size();
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
    os << "Number of nodes: " << numNodes() << "\n";
    os << "Number of edges: " << numEdges() << "\n";
    for (auto it = mConnectivity.begin(), itEnd = mConnectivity.end(); it != itEnd; ++it) {
        // Aliases for convenience
        const int&           node        = it->first;
        const std::set<int>& connections = it->second;

        // // Reverse lookup to see if this node is connected with any nodes
        // // whose index are less than it.
        // std::set<int> tmp;
        // for (auto j = mConnectivity.begin(); j != it; ++j) {
        //     // Alias
        //     const auto& s = j->second;
        //     if (s.find(node) != s.end()) {
        //         tmp.insert(j->first);
        //     }
        // }

        // // Merge what's found with the set of saved connections.
        // std::set<int> allConnections;
        // std::set_union(tmp.begin(), tmp.end(),
        //                connections.begin(), connections.end(),
        //                std::inserter(allConnections, allConnections.begin()));

        // Print
        os << "  Node: " << node;
        if (!connections.empty()) {
            os << " ---> ";
            for (int n : connections) {
                os << n << " ";
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

class FundamentalCyclesOp
{
public:

    FundamentalCyclesOp(const UndirectedGraph& graph) : mGraph(graph)
    {
        // empty
    }

    void compute()
    {
        using MapConstIterator = std::map<int, TreeNode>::const_iterator;
        using MapIterator      = std::map<int, TreeNode>::iterator;

        const std::set<int> graphNodes = mGraph.getNodes();

        // Spanning tree managed by a map because the indices of the nodes
        // in the graph need not be contiguous.
        std::map<int, TreeNode> spTree;
        for (const int n : graphNodes) {
            spTree.insert({n, TreeNode(n)});
        }

        std::stack<MapConstIterator> treeIterStack;
        treeIterStack.push(spTree.cbegin());

        // Copy the UndirectedGraph since we will be removing edges
        UndirectedGraph spGraph = mGraph;

        while (!treeIterStack.empty()) {
            const auto iter = treeIterStack.top();
            treeIterStack.pop();

            const int       currentNodeIdx = iter->first;
            const TreeNode& currentNode    = iter->second;

            // std::cout << "DEBUG ::: Processing node " << currentNodeIdx << std::endl;

            const std::set<int> connections = spGraph.getConnections(currentNodeIdx);

            for (const int j : connections) {

                // std::cout << "DEBUG :::   Other node " << j << std::endl;

                // Check if the other node is already in the spanning tree
                // by checking if its parent is null.
                const auto jIter = spTree.find(j);
                TreeNode& otherNode = jIter->second;
                if (otherNode.getParent()) {
                    // Get unique paths between both nodes
                    UndirectedGraph ga, gb;
                    _findPathToRoot(currentNode, ga);
                    _findPathToRoot(otherNode, gb);

                    // Also need to add the edge between currentNode and otherNode,
                    // but only to one of the graphs (doesn't matter which one)
                    ga.addEdge(currentNodeIdx, j);

                    // Perform symmetric difference to obtain the fundamental cycle.
                    mFundamentalCycles.push_back(ga ^ gb);

#ifdef _DEBUG_
                    std::cout << "######################################" << std::endl;
                    std::cout << "Found a cycle" << std::endl;
                    std::cout << ga << std::endl;
                    std::cout << gb << std::endl;
                    std::cout << mFundamentalCycles.back() << std::endl;
                    std::cout << "######################################" << std::endl;
#endif

                    // TODO Do I really need to do this much work just to get the
                    // fundamental cycle?
                } else {
                    // Other node is not contained in the tree, so we add it
                    // by setting its correct parent.
                    otherNode.setParent(&currentNode);
                    // Add node to the stack to be processed.
                    treeIterStack.push(jIter);
                }

                // Either way, remove this edge.
                spGraph.removeEdge(currentNodeIdx, j);

                // std::cout << spGraph << std::endl;
            }
        }
    }

    const std::vector<UndirectedGraph>& getFundamentalCycles() const { return mFundamentalCycles; }

private:

    class TreeNode
    {
    public:
        TreeNode(const int n) : mIndex(n), mParent(nullptr) {}

        int             getIndex()  const { return mIndex; }
        const TreeNode* getParent() const { return mParent; }

        void            setParent(const TreeNode* ptr) { mParent = ptr; }
    private:
        const int       mIndex;
        const TreeNode* mParent;
    };

    void _findPathToRoot(const TreeNode& node, UndirectedGraph& graph)
    {
        if (node.getParent()) {
            graph.addEdge(node.getIndex(), node.getParent()->getIndex());
            _findPathToRoot(*(node.getParent()), graph);
        }
    };

    const UndirectedGraph& mGraph;
    std::vector<UndirectedGraph> mFundamentalCycles;
};

} // namespace vspc

#endif // UNDIRECTED_GRAPH_HAS_BEEN_INCLUDED
