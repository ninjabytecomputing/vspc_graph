#ifndef SSSR_HAS_BEEN_INCLUDED
#define SSSR_HAS_BEEN_INCLUDED

#include "UndirectedGraph.h"
#include "UpperTriangularMatrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator> // std::distance
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace vspc
{

// NOTE I assume that the graphs we're working with only have one connected
// component. It might be a smart thing to implement a method that counts
// the number of components...

class SSSR
{
public:
    using NodeType = UndirectedGraph::NodeType;

    class Path
    {
    public:
        using iterator       = std::list<NodeType>::iterator;
        using const_iterator = std::list<NodeType>::const_iterator;

        Path() = default;

        Path(const NodeType a, const NodeType b) : mList({a, b}) {}

        // Default copy, move, destructor are fine.

        Path& append(const Path& path);

        iterator       begin()        { return mList.begin(); }
        const_iterator begin()  const { return mList.begin(); }
        const_iterator cbegin() const { return mList.cbegin(); }

        iterator       end()        { return mList.end(); }
        const_iterator end()  const { return mList.end(); }
        const_iterator cend() const { return mList.cend(); }

        size_t length() const { return mList.size() - 1; }

        std::string str() const;
        std::string strAsCSV() const;

    private:
        std::list<NodeType> mList;
    };

    class CandidateRing
    {
    public:
        CandidateRing(const size_t length,
                     const std::list<Path>& p)
            : mLength(length), mPath(&p), mPathP(nullptr) {}

        CandidateRing(const size_t length,
                     const std::list<Path>& p,
                     const std::list<Path>& pp)
            : mLength(length), mPath(&p), mPathP(&pp) {}

        size_t length() const { return mLength; }
        std::list<Path> const * const getPtrPaths() const { return mPath; }
        std::list<Path> const * const getPtrPathPs() const { return mPathP; }

        friend bool operator<(const CandidateRing& lhs, const CandidateRing& rhs);

    private:
        const size_t mLength;
        std::list<Path> const * const mPath;
        std::list<Path> const * const mPathP;
    };

    explicit SSSR(const UndirectedGraph& graph);

    const std::vector<UndirectedGraph>& run();

    const std::vector<UndirectedGraph>& getSSSR() const { return mSSSR; }

private:
    void _initializePID();
    void _makeCandidateRing();
    void _constructSSSR();
    void _convertNodeIndices();

    /// @{
    /// @details These are methods used for primarily used for debugging, e.g.
    /// printing out the true indices of paths/graphs that the SSSR algorithm
    /// found. Before the SSSR algorithm is complete, the indices that are used
    /// by the internal data structures won't match the indices of the input
    /// graph if (1) the indexing did not start with 0, and (2) the indices
    /// are not contiguous. To ensure the most compact data structures for this
    /// algorithm, we first perform a "transformation" on the indices so they
    /// start with 0 and are contiguous. Then, the final step of the algorithm
    /// will convert the indexing back to that of the input graph.
    void            _replaceIndices(Path& p) const;
    Path            _convertIndices(Path p) const;
    void            _replaceIndices(UndirectedGraph& g) const;
    UndirectedGraph _convertIndices(UndirectedGraph g) const;
    /// @}

    void _process(const size_t ij, const size_t ik, const size_t kj);

    bool                                   mNeedConversion;
    const UndirectedGraph&                 mGraph;
    std::vector<NodeType>                  mNodeMap;
    UpperTriangularMatrix<float>           mDold, mDnew;
    UpperTriangularMatrix<std::list<Path>> mP, mPp;
    std::list<CandidateRing>               mCandidates;
    std::vector<UndirectedGraph>           mSSSR;
};


vspc::SSSR::Path convertGraphToPath(vspc::UndirectedGraph graph)
{
    using NodeType = vspc::UndirectedGraph::NodeType;
    using Path     = vspc::SSSR::Path;

    const NodeType startNode = graph.minNode();
    NodeType node = startNode;
    std::set<NodeType> conn = graph.getConnections(node);

    Path path(node, *conn.cbegin());
    graph.removeNode(node);
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

    // Append the last path to make it a full cycle
    path.append(Path{*(--path.cend()), startNode});
    return path;
}

// -----------------------------------------------------------------------------

SSSR::Path&
SSSR::Path::append(const Path& path)
{
    if (mList.empty()) {
        mList = path.mList;
        return *this;
    }

    if (mList.back() == path.mList.front()) {
        auto it = ++(path.mList.cbegin());
        const auto itEnd = path.mList.cend();
        for (; it != itEnd; ++it) {
            mList.push_back(*it);
        }
    } else if (mList.back() == path.mList.back()) {
        auto it = ++(path.mList.crbegin());
        const auto itBegin = path.mList.crend();
        for (; it != itBegin; ++it) {
            mList.push_back(*it);
        }
    } else if (mList.front() == path.mList.front()) {
        auto it = ++(path.mList.cbegin());
        const auto itEnd = path.mList.cend();
        for (; it != itEnd; ++it) {
            mList.push_front(*it);
        }
    } else {
        auto it = ++(path.mList.crbegin());
        const auto itBegin = path.mList.crend();
        for (; it != itBegin; ++it) {
            mList.push_front(*it);
        }
    }
    return *this;
}

std::string
SSSR::Path::str() const
{
    std::ostringstream os;
    auto it = mList.cbegin();
    const auto itEnd = --(mList.cend());
    for (; it != itEnd; ++it) {
        os << *it << " -> ";
    }
    os << *it;
    return os.str();
}

std::string
SSSR::Path::strAsCSV() const
{
    std::ostringstream os;
    auto it = mList.cbegin();
    const auto itEnd = --(mList.cend());
    for (; it != itEnd; ++it) {
        os << *it << ",";
    }
    os << *it;
    return os.str();
}

SSSR::Path
merge(SSSR::Path p1, const SSSR::Path& p2)
{
    return p1.append(p2);
}

UndirectedGraph
constructGraph(const SSSR::Path& path)
{
    UndirectedGraph g;
    auto itSlow = path.cbegin();
    auto itFast = ++(path.cbegin());
    const auto itEnd = path.cend();
    while (itFast != itEnd) {
        g.addEdge(*(itSlow++), *(itFast++));
    }
    return g;
}

std::ostream&
operator<<(std::ostream& os, const SSSR::Path& obj)
{
    return os << obj.str();
}

// -----------------------------------------------------------------------------

bool operator<(const SSSR::CandidateRing& lhs, const SSSR::CandidateRing& rhs)
{
    return lhs.mLength < rhs.mLength;
}

// -----------------------------------------------------------------------------

SSSR::SSSR(const UndirectedGraph& graph)
    : mNeedConversion(false)
    , mGraph(graph)
    , mNodeMap(graph.numNodes())
    , mDold(graph.numNodes(), INFINITY)
    , mDnew(graph.numNodes(), INFINITY)
    , mP(graph.numNodes())
    , mPp(graph.numNodes())
{
    const std::set<NodeType> nodes = graph.getNodes();
    std::copy(nodes.cbegin(), nodes.cend(), mNodeMap.begin());

    // Scan node indices to see a conversion is required later
    bool contiguous = true;
    for (size_t i = 0, n = nodes.size() - 1; i < n; ++i) {
        if (mNodeMap[i+1] - mNodeMap[i] != 1) {
            contiguous = false;
            break;
        }
    }
    if (mNodeMap[0] != 0 || !contiguous) {
        mNeedConversion = true;
    }

    for (size_t i = 0, n = nodes.size(); i < n; ++i) {
        for (NodeType jj : graph.getConnections(mNodeMap[i])) {
            auto it = std::lower_bound(mNodeMap.cbegin(), mNodeMap.cend(), jj);
            const size_t j = std::distance(mNodeMap.cbegin(), it);
            const size_t ij = mDold.index(i, j);
            mDold(ij) = mDnew(ij) = 1.f;
            mP(ij) = {Path(i, j)};
        }
    }
}

const std::vector<UndirectedGraph>&
SSSR::run()
{
    _initializePID();
    _makeCandidateRing();
    _constructSSSR();
    _convertNodeIndices();
    return mSSSR;
}

void
SSSR::_initializePID()
{
    for (size_t k = 0, n = mNodeMap.size(); k < n; ++k) {

        // i < k
        for (size_t i = 0; i < k; ++i) {
            const size_t ik = mDold.index(i, k);

            // i < j < k
            for (size_t j = i + 1; j < k; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(j, k);

                _process(ij, ik, kj);
            }
        }

        // i < k
        for (size_t i = 0; i < k; ++i) {
            const size_t ik = mDold.index(i, k);

            // k <= j
            for (size_t j = k; j < n; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(k, j);

                _process(ij, ik, kj);
            }
        }

        // k <= i
        for (size_t i = k; i < n; ++i) {
            const size_t ik = mDold.index(k, i);

            // i < j
            for (size_t j = i + 1; j < n; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(k, j);

                _process(ij, ik, kj);
            }
        }

        mDold = mDnew;
    }

#ifdef _VERBOSE
    std::cout << "Matrix P" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            std::cout << "  (" << i << ", " << j << ") ::: " << mP(i, j).size() << "\n";
            for (const auto& path : mP(i, j)) {
                std::cout << "    " << path << "\n";
            }
        }
        std::cout << "\n";
    }
    std::cout << "Matrix P'" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            std::cout << "  (" << i << ", " << j << ") ::: " << mPp(i, j).size() << "\n";
            for (const auto& path : mPp(i, j)) {
                std::cout << "    " << path << "\n";
            }
        }
        std::cout << "\n";
    }
    std::cout << "Distance" << std::endl;
    std::cout << mDnew << std::endl;
#endif
}

void
SSSR::_makeCandidateRing()
{
    for (size_t i = 0, n = mNodeMap.size(); i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const size_t ij = mDnew.index(i, j);

            if (std::isinf(mDnew(ij)) || (mP(ij).size() == 1 && mPp(ij).empty())) {
                continue;
            } else {
                if (!mPp(ij).empty()) {
                    mCandidates.emplace_back(2 * mDnew(ij) + 1, mP(ij), mPp(ij));
                } else {
                    mCandidates.emplace_back(2 * mDnew(ij), mP(ij));
                }
            }
        }
    }
    // Sort so we construct our SSSR set starting with the shorted cycles.
    mCandidates.sort();
}

void
SSSR::_constructSSSR()
{
    // Functor used to check whether the newly generated graph is a cycle that
    // should be inserted in our SSSR, returning true if so.
    auto checkCycle = [this](const UndirectedGraph& g,
                             const Path& p1,
                             const Path& p2)
    {
        // Edges overlapped, not a good candidate
        if (g.numEdges() != p1.length() + p2.length()) { return false; }
        // Nodes overlapped, not a good candidate
        if (g.numNodes() != g.numEdges()) { return false; }

        // Check if we've already generated this cycle by simply checking
        // for equality against what's already generated.
        for (const UndirectedGraph& ring : mSSSR) {
            if (g == ring) { return false; }
        }
        return true;
    };

    for (const CandidateRing& candidate : mCandidates) {
        if (candidate.length() % 2 == 1) {
#ifdef _DEBUG
            assert(candidate.getPtrPathPs());
#endif
            for (const Path& lp : *candidate.getPtrPathPs()) {
                for (const Path& sp : *candidate.getPtrPaths()) {

                    UndirectedGraph g0 = constructGraph(lp);
                    UndirectedGraph g1 = constructGraph(sp);
                    g0 ^= g1;

                    if (checkCycle(g0, lp, sp)) {
                        mSSSR.push_back(g0);
                    }
                }
            }
        } else {
            const std::list<Path>& paths = *candidate.getPtrPaths();
            const auto itEnd = paths.cend();
            for (auto itA = paths.cbegin(); itA != itEnd; ++itA) {
                auto itBInit = itA;
                for (auto itB = ++itBInit; itB != itEnd; ++itB) {
                    const Path& p = *itA;
                    const Path& q = *itB;

                    UndirectedGraph g0 = constructGraph(p);
                    UndirectedGraph g1 = constructGraph(q);
                    g0 ^= g1;

                    if (checkCycle(g0, p, q)) {
                        mSSSR.push_back(g0);
                    }
                }
            }
        }
    }
}

void
SSSR::_convertNodeIndices()
{
    // Check for early out
    if (!mNeedConversion) return;

    std::vector<UndirectedGraph> result(mSSSR.size());
    for (size_t i = 0, n = mSSSR.size(); i < n; ++i) {
        for (const auto& [node, conn] : mSSSR[i]) {
            for (const NodeType c : conn) {
                result[i].addEdge(mNodeMap[node], mNodeMap[c]);
            }
        }
    }
    mSSSR = result;
}

void
SSSR::_replaceIndices(Path& p) const
{
    for (NodeType& n : p) {
        n = mNodeMap[n];
    }
}

SSSR::Path
SSSR::_convertIndices(Path p) const
{
    _replaceIndices(p);
    return p;
}

void
SSSR::_replaceIndices(UndirectedGraph& g) const
{
    UndirectedGraph tmp;
    for (const auto& [node, conn] : g) {
        for (const NodeType c : conn) {
            tmp.addEdge(mNodeMap[node], mNodeMap[c]);
        }
    }
    g = tmp;
}

UndirectedGraph
SSSR::_convertIndices(UndirectedGraph g) const
{
    _replaceIndices(g);
    return g;
}

void
SSSR::_process(const size_t ij, const size_t ik, const size_t kj)
{
    if (mDold(ij) > mDold(ik) + mDold(kj)) {
        if (!std::isinf(mDold(ij)) && mDold(ij) == mDold(ik) + mDold(kj) + 1) {
#ifdef _DEBUG
            assert(!mP(ij).empty());
#endif
            mPp(ij) = {mP(ij).front()};
        } else {
            mPp(ij).clear();
        }
        mDnew(ij) = mDold(ik) + mDold(kj);
#ifdef _DEBUG
        assert(!mP(ik).empty());
        assert(!mP(kj).empty());
#endif
        mP(ij) = {merge(mP(ik).front(), mP(kj).front())};
    } else if (!std::isinf(mDold(ij)) && mDold(ij) == mDold(ik) + mDold(kj)) {
#ifdef _DEBUG
        assert(!mP(ik).empty());
        assert(!mP(kj).empty());
#endif
        mP(ij).push_back(merge(mP(ik).front(), mP(kj).front()));
    } else if (!std::isinf(mDold(ij)) && mDold(ij) == mDold(ik) + mDold(kj) - 1) {
#ifdef _DEBUG
        assert(!mP(ik).empty());
        assert(!mP(kj).empty());
#endif
        mPp(ij).push_back(merge(mP(ik).front(), mP(kj).front()));
    } else {
        mDnew(ij) = mDold(ij);
    }
}

} // namespace vspc

#endif // SSSR_HAS_BEEN_INCLUDED
