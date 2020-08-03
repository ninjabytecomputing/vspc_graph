#ifndef SSSR_HAS_BEEN_INCLUDED
#define SSSR_HAS_BEEN_INCLUDED

#include "UndirectedGraph.h"
#include "UpperTriangularMatrix.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <iterator> // std::distance
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/blocked_range2d.h>

namespace vspc
{

// NOTE I assume that the graphs we're working with only have one connected
// component. It might be a smart thing to implement a method that counts
// the number of components...

class SSSR
{
public:
    using NodeType = UndirectedGraph::NodeType;

    class Settings
    {
    public:
        Settings() {};
        void setMaxCycleLength(const int len) { mMaxCycleLength = len; }
        int getMaxCycleLength() const { return mMaxCycleLength; }
    private:
        int mMaxCycleLength = -1;
    };

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

    explicit SSSR(const UndirectedGraph& graph,
                  const Settings& settings = {});

    const std::vector<UndirectedGraph>& run();

    const std::vector<UndirectedGraph>& getSSSR() const { return mSSSR; }

private:
    void _initializePID();
    void _makeCandidateRings();
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

    UpperTriangularMatrix<float>           mDold, mDnew;
    UpperTriangularMatrix<std::list<Path>> mP, mPp;
    std::list<CandidateRing>               mCandidates;
    std::vector<UndirectedGraph>           mSSSR;
    std::vector<NodeType>                  mNodeMap;
    const UndirectedGraph&                 mGraph;
    const Settings&                        mSettings;
    bool                                   mNeedConversion;
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

SSSR::SSSR(const UndirectedGraph& graph, const Settings& settings)
    : mDold(graph.numNodes(), INFINITY)
    , mDnew(graph.numNodes(), INFINITY)
    , mP(graph.numNodes())
    , mPp(graph.numNodes())
    , mNodeMap(graph.numNodes())
    , mGraph(graph)
    , mSettings(settings)
    , mNeedConversion(false)
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
    std::cout << "Initialize" << std::endl;
    _initializePID();

    std::cout << "Candidates" << std::endl;
    _makeCandidateRings();

    std::cout << "Construct" << std::endl;
    _constructSSSR();

    std::cout << "Convert" << std::endl;
    _convertNodeIndices();

    return mSSSR;
}

void
SSSR::_initializePID()
{
    const size_t numElem = mDold.numElements();
    const size_t dim     = mDold.dim();

    // Functor for computing the offset of the first dimension of indices
    auto firstDimOffset = [numElem, dim](size_t i) {
        return numElem - ((dim - i) * (dim + 1 - i) / 2 + i);
    };

    for (size_t k = 0, n = mNodeMap.size(); k < n; ++k) {
        const size_t kOffset = firstDimOffset(k);

        auto t1 = std::chrono::high_resolution_clock::now();

        tbb::parallel_invoke(
            [this, &firstDimOffset, k] {
#if 1
                // i < k
                tbb::parallel_for(tbb::blocked_range<size_t>(0, k),
                    [this, &firstDimOffset, k](tbb::blocked_range<size_t>& r)
                    {
                        for (auto i = r.begin(), iEnd = r.end(); i < iEnd; ++i) {
                            const size_t iOffset = firstDimOffset(i);
                            const size_t ik      = iOffset + k;

                            // i < j < k
                            for (size_t j = i + 1; j < k; ++j) {
                                // const size_t ij = mDold.index(i, j);
                                // const size_t kj = mDold.index(j, k);
                                const size_t ij = iOffset + j;
                                const size_t kj = firstDimOffset(j) + k;
                                this->_process(ij, ik, kj);
                            }
                        }
                    }
                );
#else
                // i < k
                for (size_t i = 0; i < k; ++i) {
                    // const size_t ik = mDold.index(i, k);
                    const size_t iOffset = firstDimOffset(i);
                    const size_t ik      = iOffset + k;

                    // i < j < k
                    for (size_t j = i + 1; j < k; ++j) {
                        // const size_t ij = mDold.index(i, j);
                        // const size_t kj = mDold.index(j, k);
                        const size_t ij = iOffset + j;
                        const size_t kj = firstDimOffset(j) + k;
                        this->_process(ij, ik, kj);
                    }
                }
#endif
            },
            [this, &firstDimOffset, k, kOffset, n] {
#if 1
                // i < k, k < j
                tbb::parallel_for(tbb::blocked_range2d<size_t, size_t>(0, k, k + 1, n),
                    [this, &firstDimOffset, k, kOffset](tbb::blocked_range2d<size_t, size_t>& r)
                    {
                        auto rows = r.rows();
                        auto cols = r.cols();
                        size_t iOffset, ij, ik, kj;
                        for (auto i = rows.begin(), iEnd = rows.end(); i < iEnd; ++i) {
                            iOffset = firstDimOffset(i);
                            ik      = iOffset + k;
                            for (auto j = cols.begin(), jEnd = cols.end(); j < jEnd; ++j) {
                                ij = iOffset + j;
                                kj = kOffset + j;
                                this->_process(ij, ik, kj);
                            }
                        }
                    }
                );
#else
                // i < k
                for (size_t i = 0; i < k; ++i) {
                    // const size_t ik = mDold.index(i, k);
                    const size_t iOffset = firstDimOffset(i);
                    const size_t ik      = iOffset + k;

                    // k < j
                    for (size_t j = k + 1; j < n; ++j) {
                        // const size_t ij = mDold.index(i, j);
                        // const size_t kj = mDold.index(k, j);
                        const size_t ij = iOffset + j;
                        const size_t kj = kOffset + j;
                        _process(ij, ik, kj);
                    }
                }
#endif
            },
            [this, &firstDimOffset, k, kOffset, n] {
#if 1
                // k < i
                tbb::parallel_for(tbb::blocked_range<size_t>(k + 1, n),
                    [this, &firstDimOffset, kOffset, n](tbb::blocked_range<size_t>& r)
                    {
                        for (auto i = r.begin(), iEnd = r.end(); i < iEnd; ++i) {
                            const size_t iOffset = firstDimOffset(i);
                            const size_t ik      = kOffset + i;

                            // i < j
                            for (size_t j = i + 1; j < n; ++j) {
                                // const size_t ij = mDold.index(i, j);
                                // const size_t kj = mDold.index(k, j);
                                const size_t ij = iOffset + j;
                                const size_t kj = kOffset + j;
                                this->_process(ij, ik, kj);
                            }
                        }
                    }
                );
#else
                // k < i
                for (size_t i = k + 1; i < n; ++i) {
                    // const size_t ik = mDold.index(k, i);
                    const size_t iOffset = firstDimOffset(i);
                    const size_t ik      = kOffset + i;

                    // i < j
                    for (size_t j = i + 1; j < n; ++j) {
                        // const size_t ij = mDold.index(i, j);
                        // const size_t kj = mDold.index(k, j);
                        const size_t ij = iOffset + j;
                        const size_t kj = kOffset + j;
                        _process(ij, ik, kj);
                    }
                }
#endif
            }
        );

        mDold = mDnew;

        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        std::cout << k << " ::: " << duration << std::endl;
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
SSSR::_makeCandidateRings()
{
    for (size_t i = 0, n = mNodeMap.size(); i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const size_t ij = mDnew.index(i, j);

            if (std::isinf(mDnew(ij)) || (mP(ij).size() == 1 && mPp(ij).empty())) {
                continue;
            } else {
                if (mP(ij).size() > 1) {
                    const size_t len = 2 * mDnew(ij);
                    if (len <= mSettings.getMaxCycleLength()) {
                        mCandidates.emplace_back(len, mP(ij));
                    }
                } else {
                    const size_t len = 2 * mDnew(ij) + 1;
                    if (len <= mSettings.getMaxCycleLength()) {
                        mCandidates.emplace_back(len, mP(ij), mPp(ij));
                    }
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

        // Check if this new cycle satisfies one of the following:
        //   1) already generated, which only needs to be checked against
        //      cycles of the same length; or
        //   2) contains all of a cycle that has already been generated.
        // If either are true, we return false, i.e. we don't insert this
        // new cycle into the SSSR set.
        const std::set<NodeType> gNodes = g.getNodes();
        for (const UndirectedGraph& ring : mSSSR) {
            if (ring.numEdges() == g.numEdges()) {
                if (g == ring) { return false; }
            } else {
                const std::set<NodeType> rNodes = ring.getNodes();
                if (std::includes(gNodes.cbegin(), gNodes.cend(),
                                  rNodes.cbegin(), rNodes.cend())) { return false; }
            }

        }
        return true;
    };

    for (const CandidateRing& candidate : mCandidates) {
        if (candidate.length() % 2 == 1) {
#ifdef _DEBUG
            assert(candidate.getPtrPathPs());
            assert(candidate.getPtrPaths());
            assert(candidate.getPtrPaths()->size() == 1);
#endif
            // Only one short path
            const Path& sp           = candidate.getPtrPaths()->front();
            const UndirectedGraph g1 = constructGraph(sp);

            for (const Path& lp : *candidate.getPtrPathPs()) {
                UndirectedGraph g = constructGraph(lp);
                g ^= g1;

                if (checkCycle(g, lp, sp)) {
                    mSSSR.push_back(g);
                }
            }
        } else {
#ifdef _DEBUG
            assert(candidate.getPtrPaths());
#endif
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
    float dSum = mDold(ik) + mDold(kj);
    if (dSum > mSettings.getMaxCycleLength() / 2) {
        dSum = INFINITY;
    }
    const bool isNotInf = !std::isinf(mDold(ij));

    if (mDold(ij) > dSum) {
        if (isNotInf && mDold(ij) == dSum + 1) {
#ifdef _DEBUG
            assert(!mP(ij).empty());
#endif
            mPp(ij) = {mP(ij).front()};
        } else {
            mPp(ij).clear();
        }
        mDnew(ij) = dSum;
#ifdef _DEBUG
        assert(!mP(ik).empty());
        assert(!mP(kj).empty());
#endif
        mP(ij) = {merge(mP(ik).front(), mP(kj).front())};
    } else if (isNotInf && mDold(ij) == dSum) {
#ifdef _DEBUG
        assert(!mP(ik).empty());
        assert(!mP(kj).empty());
#endif
        mP(ij).push_back(merge(mP(ik).front(), mP(kj).front()));
    } else if (isNotInf && mDold(ij) == dSum - 1) {
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
