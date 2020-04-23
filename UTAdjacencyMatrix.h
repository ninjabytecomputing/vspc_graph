#ifndef UTADJACENCY_MATRIX_HAS_BEEN_INCLUDED
#define UTADJACENCY_MATRIX_HAS_BEEN_INCLUDED

#include <cassert>
#include <string>
#include <vector>

namespace vspc
{

class UTAdjacencyMatrix
{
public:
    explicit UTAdjacencyMatrix(const size_t numNodes);

    // Default copy, move, and destructor are fine.

    void addEdge(const size_t i, const size_t j);
    void removeEdge(const size_t i, const size_t j);

    bool operator()(const size_t i, const size_t j);

    size_t getnNumNodes() const { return mNumNodes; }
    size_t getNumEdges()  const { return mNumEdges; }

    UTAdjacencyMatrix& operator^=(const UTAdjacencyMatrix& rhs);

    // TODO
    std::string str() const;

private:
    // TODO
    size_t _index(const size_t i, const size_t j) const;

    std::vector<bool> mElems;
    const size_t      mNumNodes;
    size_t            mNumEdges;
    const size_t      mIndexFactor;
};

UTAdjacencyMatrix operator^(UTAdjacencyMatrix lhs,
                            const UTAdjacencyMatrix& rhs);

// -----------------------------------------------------------------------------

UTAdjacencyMatrix::UTAdjacencyMatrix(const size_t numNodes)
    : mElems((numNodes * (numNodes - 1)) / 2)
    , mNumNodes(numNodes)
    , mNumEdges(0)
    , mIndexFactor(1 + 2 * (numNodes - 2))
{
    // empty
}

void
UTAdjacencyMatrix::addEdge(const size_t i, const size_t j)
{
    if (const auto idx = _index(i, j); !mElems[idx]) {
        mElems[idx] = true;
        ++mNumEdges;
    }
}

bool
UTAdjacencyMatrix::operator()(const size_t i, const size_t j)
{
    return (i == j ? true : mElems[_index(i, j)]);
}

void
UTAdjacencyMatrix::removeEdge(const size_t i, const size_t j)
{
    if (const auto idx = _index(i, j); mElems[idx]) {
        mElems[idx] = false;
        --mNumEdges;
    }
}

UTAdjacencyMatrix&
UTAdjacencyMatrix::operator^=(const UTAdjacencyMatrix& rhs)
{
    // TODO Wrap assert in debug
    assert(mNumNodes == rhs.mNumNodes);

    mNumEdges = 0;
    for (size_t i = 0, n = mElems.size(); i < n; ++i) {
        if (mElems[i] != rhs.mElems[i]) {
            mElems[i] = true;
            +mNumEdges;
        } else {
            mElems[i] = false;
        }
    }
    return *this;
}

size_t
UTAdjacencyMatrix::_index(const size_t i, const size_t j) const
{
    // TODO Wrap asserts in debug
    assert(i < mNumNodes);
    assert(j < mNumNodes);
    assert(i != j);

    return i;
}

UTAdjacencyMatrix operator^(UTAdjacencyMatrix lhs,
                            const UTAdjacencyMatrix& rhs)
{
    return lhs ^= rhs;
}

} // namespace vspc

#endif // UTADJACENCY_MATRIX_HAS_BEEN_INCLUDED
