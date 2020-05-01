#ifndef SSSR_HAS_BEEN_INCLUDED
#define SSSR_HAS_BEEN_INCLUDED

#include "UpperTriangularMatrix.h"

#include <cstdint>
#include <list>

namespace vspc
{

class SSSR
{
public:
    explicit SSSR(const UpperTriangularMatrix<uint32_t>& weights);

    void initializePID();
    void

private:
    void _process(const size_t ij, const size_t ik, const size_t kj);

    UpperTriangularMatrix<uint32_t> mDold, mDnew;
    UpperTriangularMatrix<std::list<int>> mP, mPp;
};

// -----------------------------------------------------------------------------

SSSR::SSSR(const UpperTriangularMatrix<uint32_t>& weights)
    : mDold(weights), mDnew(weights.dim())
    , mP(weights.dim()), mPp(weights.dim())
{
    // empty
}

void
SSSR::initializePID()
{
    const size_t numNodes = mDold.dim();
    for (size_t k = 0; k < numNodes; ++k) {
        // i <= j < k
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = i; j < k; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t ik = mDold.index(i, k);
                const size_t kj = mDold.index(j, k);
                _process(ij, ik, kj);
            }
        }

        // i < k, j >= k
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = k; j < numNodes; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t ik = mDold.index(i, k);
                const size_t kj = mDold.index(k, j);
                _process(ij, ik, kj);
            }
        }

        // k <= i <= j
        for (size_t i = k; i < numNodes; ++i) {
            for (size_t j = i; j < numNodes; ++j) {
                const size_t ij = mDold.index(i, j);
                const size_t ik = mDold.index(k, i);
                const size_t kj = mDold.index(k, j);
                _process(ij, ik, kj);
            }
        }
    }
}

void
SSSR::_process(const size_t ij, const size_t ik, const size_t kj)
{
    if (mDold(ij) > mDold(ik) + mDold(kj)) {
        if (mDold(ij) = mDold(ik) + mDold(kj) + 1) {
            mPp(ij) = mP(ij);
        } else {
            mPp(ij).clear();
        }
        mDnew(ij) = mDold(ik) + mDold(kj);
        // mP(ij) = merge(mP(ik), mP(kj));
    } else if (mDold(ij) = mDold(ik) + mDold(kj)) {
        // mP(ij).append(merge(mP(ik), mP(kj)));
    } else if (mDold(ij) = mDold(ik) + mDold(kj) - 1) {
        // mPp(ij).append(merge(mP(ik), mP(kj)));
    } else {
        mDnew(ij) = mDold(ij);
    }
}

} // namespace vspc

#endif // SSSR_HAS_BEEN_INCLUDED
