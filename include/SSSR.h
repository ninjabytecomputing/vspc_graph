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

    void _process(const size_t ij, const size_t ik, const size_t kj);

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
    : mGraph(graph)
    , mNodeMap(graph.numNodes())
    , mDold(graph.numNodes(), INFINITY)
    , mDnew(graph.numNodes(), INFINITY)
    , mP(graph.numNodes())
    , mPp(graph.numNodes())
{
    const std::set<NodeType> nodes = graph.getNodes();
    std::copy(nodes.cbegin(), nodes.cend(), mNodeMap.begin());

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
    for (const CandidateRing& candidate : mCandidates) {
        if (candidate.length() % 2 == 1) {
#ifdef _DEBUG
            assert(candidate.getPtrPathPs());
#endif
            for (const Path& longPath : *candidate.getPtrPathPs()) {
                const Path& shortPath = candidate.getPtrPaths()->front();
                for (const Path& shortPath : *candidate.getPtrPaths()) {

                    UndirectedGraph g0 = constructGraph(longPath);
                    UndirectedGraph g1 = constructGraph(shortPath);
                    g0 ^= g1;

                    // Edges overlapped, not a good candidate
                    if (g0.numEdges() != longPath.length() + shortPath.length()) {
                        continue;
                    }

                    // Nodes overlapped, not a good candidate
                    if (g0.numNodes() != g0.numEdges()) {
                        continue;
                    }

                    NodeType a = *shortPath.cbegin();
                    NodeType b = *(--shortPath.cend());
                    bool alreadyConnected = false;
                    if (a > b) std::swap(a, b);
                    for (const UndirectedGraph& ring : mSSSR) {
                        if (ring.hasNode(a) && ring.hasNode(b)) {
                            alreadyConnected = true;
                            break;
                        }
                    }

                    if (!alreadyConnected) {
                        mSSSR.push_back(g0);
                    }
                }
            }
        } else {
            const std::list<Path>& shortPaths = *candidate.getPtrPaths();
            const auto itEnd = shortPaths.cend();
            for (auto itA = shortPaths.cbegin(); itA != itEnd; ++itA) {
                auto itBInit = itA;
                for (auto itB = ++itBInit; itB != itEnd; ++itB) {
                    const Path& p1 = *itA;
                    const Path& p2 = *itB;

                    UndirectedGraph g0 = constructGraph(p1);
                    UndirectedGraph g1 = constructGraph(p2);
                    g0 ^= g1;

                    // Edges overlapped, not a good candidate
                    if (g0.numEdges() != p1.length() + p2.length()) {
                        continue;
                    }

                    // Nodes overlapped, not a good candidate
                    if (g0.numNodes() != g0.numEdges()) {
                        continue;
                    }

                    NodeType a = *p1.cbegin();
                    NodeType b = *(--p1.cend());
                    bool alreadyConnected = false;
                    if (a > b) std::swap(a, b);
                    for (const UndirectedGraph& ring : mSSSR) {
                        if (ring.hasNode(a) && ring.hasNode(b)) {
                            alreadyConnected = true;
                            break;
                        }
                    }

                    if (!alreadyConnected) {
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
