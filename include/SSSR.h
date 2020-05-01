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

namespace vspc
{

class Path
{
public:
    Path() = default;

    Path(const size_t a, const size_t b)
    {
        mList.push_back(a);
        mList.push_back(b);
    }

    // Default copy, move, destructor are fine.

    Path& append(const Path& path)
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

    size_t length() const { return mList.size(); }

    const std::list<size_t>& nodes() const { return mList; }

    std::string str() const
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

    std::string rstr() const
    {
        std::ostringstream os;
        auto it = mList.crbegin();
        const auto itEnd = --(mList.crend());
        for (; it != itEnd; ++it) {
            os << *it << " -> ";
        }
        os << *it;
        return os.str();
    }

private:
    std::list<size_t> mList;
};

std::ostream&
operator<<(std::ostream& os, const Path& obj)
{
    return os << obj.str();
}

Path
merge(Path p1, const Path& p2)
{
    return p1.append(p2);
}

// -----------------------------------------------------------------------------

class SSSR
{
public:
    explicit SSSR(const UndirectedGraph& graph);

    void initializePID();
    void makeCandidateSet();

private:
    void _process(const size_t ij, const size_t ik, const size_t kj);

    UpperTriangularMatrix<float>           mDold, mDnew;
    UpperTriangularMatrix<std::list<Path>> mP, mPp;
};

// -----------------------------------------------------------------------------

SSSR::SSSR(const UndirectedGraph& graph)
    : mDold(graph.numNodes(), INFINITY)
    , mDnew(graph.numNodes(), INFINITY)
    , mP(graph.numNodes())
    , mPp(graph.numNodes())
{
    const std::set<int> nodes = graph.getNodes();
    for (auto i : nodes) {
        for (auto&& j : graph.getConnections(i)) {
            const size_t ij = mDold.index(i, j);
            mDold(ij) = 1.f;
            mDnew(ij) = 1.f;
            mP(ij) = {Path(i, j)};
        }
    }
}

void
SSSR::initializePID()
{
    const auto numNodes = mDold.dim();
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

#ifdef _DEBUG
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
SSSR::makeCandidateSet()
{
    const size_t n = mDnew.dim();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const size_t ij = mDnew.index(i, j);
            if (std::isinf(mDnew(ij)) || (mP(ij).size() == 1 && mPp(ij).empty())) {
                continue;
            } else {
                if (!mPp(ij).empty()) {
                    std::cout << "Candidate odd ring" << std::endl;
                } else {
                    std::cout << "Candidate even ring" << std::endl;
                }
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
