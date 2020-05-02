#ifndef SSSR_HAS_BEEN_INCLUDED
#define SSSR_HAS_BEEN_INCLUDED

#include "UndirectedGraph.h"
#include "UpperTriangularMatrix.h"

#include <cmath>
#include <cstdint>
#include <iostream>
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
    class Path
    {
    public:
        using iterator       = std::list<size_t>::iterator;
        using const_iterator = std::list<size_t>::const_iterator;

        Path() = default;

        Path(const size_t a, const size_t b) : mList({a, b}) {}

        // Default copy, move, destructor are fine.

        Path& append(const Path& path);

        iterator       begin()        { return mList.begin(); }
        const_iterator begin()  const { return mList.begin(); }
        const_iterator cbegin() const { return mList.cbegin(); }

        iterator       end()        { return mList.end(); }
        const_iterator end()  const { return mList.end(); }
        const_iterator cend() const { return mList.cend(); }

        size_t length() const { return mList.size(); }

        std::string str() const;

    private:
        std::list<size_t> mList;
    };

    class CandidateSet
    {
    public:
        CandidateSet(const size_t length,
                     const std::list<Path>& p,
                     const std::list<Path>& pp)
            : mLength(length), mPath(p), mPathP(pp) {}

        size_t length() const { return mLength; }
        const std::list<Path>& getPaths() const { return mPath; }
        const std::list<Path>& getPathPs() const { return mPathP; }

        friend bool operator<(const CandidateSet& lhs, const CandidateSet& rhs);

    private:
        const size_t mLength;
        const std::list<Path>& mPath;
        const std::list<Path>& mPathP;
    };

    explicit SSSR(const UndirectedGraph& graph);

    const std::vector<UndirectedGraph>& run();

    const std::vector<UndirectedGraph>& getSSSR() const { return mSSSR; }

private:
    void _initializePID();
    void _makeCandidateSet();
    void _constructSSSR();

    void _process(const size_t ij, const size_t ik, const size_t kj);

    size_t                                 mNumSSSR;
    UpperTriangularMatrix<float>           mDold, mDnew;
    UpperTriangularMatrix<std::list<Path>> mP, mPp;
    std::list<CandidateSet>                mCandidates;
    std::vector<UndirectedGraph>           mSSSR;
};

// -----------------------------------------------------------------------------

SSSR::Path&
SSSR::Path::append(const Path& path)
{
    if (mList.front() == path.mList.front()) {
        auto it = ++(path.mList.cbegin());
        const auto itEnd = path.mList.cend();
        for (; it != itEnd; ++it) {
            mList.push_front(*it);
        }
    } else if (mList.front() == path.mList.back()) {
        auto it = ++(path.mList.crbegin());
        const auto itBegin = path.mList.crend();
        for (; it != itBegin; ++it) {
            mList.push_front(*it);
        }
    } else if (mList.back() == path.mList.front()) {
        auto it = ++(path.mList.cbegin());
        const auto itEnd = path.mList.cend();
        for (; it != itEnd; ++it) {
            mList.push_back(*it);
        }
    } else {
        auto it = ++(path.mList.crbegin());
        const auto itBegin = path.mList.crend();
        for (; it != itBegin; ++it) {
            mList.push_back(*it);
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

bool operator<(const SSSR::CandidateSet& lhs, const SSSR::CandidateSet& rhs)
{
    return lhs.mLength < rhs.mLength;
}

// -----------------------------------------------------------------------------

SSSR::SSSR(const UndirectedGraph& graph)
    : mNumSSSR(graph.numEdges() - graph.numNodes() + 1)
    , mDold(graph.numNodes(), INFINITY)
    , mDnew(graph.numNodes(), INFINITY)
    , mP(graph.numNodes())
    , mPp(graph.numNodes())
{
    mSSSR.reserve(mNumSSSR);

    const std::set<int> nodes = graph.getNodes();
    for (int i : nodes) {
        for (int j : graph.getConnections(i)) {
            const size_t ij = mDold.index(i, j);
            mDold(ij) = 1.f;
            mDnew(ij) = 1.f;
            mP(ij) = {Path(i, j)};
        }
    }
}

const std::vector<UndirectedGraph>&
SSSR::run()
{
    _initializePID();
    _makeCandidateSet();
    _constructSSSR();
    return mSSSR;
}

void
SSSR::_initializePID()
{
    const size_t  numNodes = mDold.dim();
    for (size_t k = 0; k < numNodes; ++k) {
        // i <= j < k
        for (size_t i = 0; i < k; ++i) {
            const size_t ik = mDold.index(i, k);
            for (size_t j = i+1; j < k; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(j, k);
                _process(ij, ik, kj);
            }
        }

        // i < k, j >= k
        for (size_t i = 0; i < k; ++i) {
            const size_t ik = mDold.index(i, k);
            for (size_t j = k; j < numNodes; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(k, j);
                _process(ij, ik, kj);
            }
        }

        // k <= i <= j
        for (size_t i = k; i < numNodes; ++i) {
            const size_t ik = mDold.index(k, i);
            for (size_t j = i+1; j < numNodes; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t kj = mDold.index(k, j);
                _process(ij, ik, kj);
            }
        }

        mDold = mDnew;
    }

#ifdef _VERBOSE
    std::cout << "Matrix P" << std::endl;
    for (size_t i = 0; i < mP.dim(); ++i) {
        for (size_t j = i; j < mP.dim(); ++j) {
            std::cout << "  (" << i << ", " << j << ") ::: " << mP(i, j).size() << "\n";
            for (const auto& path : mP(i, j)) {
                std::cout << "    " << path << "\n";
            }
        }
        std::cout << "\n";
    }
    std::cout << "Matrix P'" << std::endl;
    for (size_t i = 0; i < mPp.dim(); ++i) {
        for (size_t j = i; j < mPp.dim(); ++j) {
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
SSSR::_makeCandidateSet()
{
    const size_t n = mDnew.dim();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const size_t ij = mDnew.index(i, j);
            if (std::isinf(mDnew(ij)) || (mP(ij).size() == 1 && mPp(ij).empty())) {
                continue;
            } else {
                if (!mPp(ij).empty()) {
                    mCandidates.emplace_back(2 * mDnew(ij) + 1, mP(ij), mPp(ij));
                } else {
                    mCandidates.emplace_back(2 * mDnew(ij), mP(ij), mPp(ij));
                }
            }
        }
    }
    mCandidates.sort();
}

void
SSSR::_constructSSSR()
{
    size_t numRings = 0;
    for (const CandidateSet& candidate : mCandidates) {
        if (candidate.length() % 2 == 1) {
            for (const Path& longPath : candidate.getPathPs()) {
                const Path& shortPath = candidate.getPaths().front();

                UndirectedGraph g = constructGraph(merge(shortPath, longPath));

                for (const UndirectedGraph& sssr : mSSSR) {
                    const auto gIt = g.cbegin();
                    const auto sssrIt = sssr.cbegin();

                    if (gIt->first == sssrIt->first && gIt->second == sssrIt->second) {
                        // The cycle we just constructed contains a cycle
                        // in our SSSR set. Perform XOR to remove the
                        // overlapping portion.
                        g ^= sssr;
                    }
                }

                // At this point, g is now cycle we haven't seen before, in
                // which case we should insert it into our SSSR set, or g
                // is empty.

                if (g.numEdges() != 0) {
                    mSSSR.push_back(g);
                    ++numRings;
                }

                if (numRings == mNumSSSR) return;
            }
        } else {
            const std::list<Path>& shortPaths = candidate.getPaths();
            auto it = shortPaths.cbegin();
            const auto itEnd = --(shortPaths.cend());

            while (it != itEnd) {
                // Why can't I just pass the iterators into merge?!
                const Path& p1 = *it;
                const Path& p2 = *(++it);

                UndirectedGraph g = constructGraph(merge(p1, p2));
                for (const UndirectedGraph& sssr : mSSSR) {
                    const auto gIt = g.cbegin();
                    const auto sssrIt = sssr.cbegin();

                    if (gIt->first == sssrIt->first && gIt->second == sssrIt->second) {
                        // The cycle we just constructed contains a cycle
                        // in our SSSR set. Perform XOR to remove the
                        // overlapping portion.
                        g ^= sssr;
                    }
                }

                // At this point, g is now cycle we haven't seen before, in
                // which case we should insert it into our SSSR set, or g
                // is empty.

                if (g.numEdges() != 0) {
                    mSSSR.push_back(g);
                    ++numRings;
                }

                if (numRings == mNumSSSR) return;
            }
        }
    }
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
