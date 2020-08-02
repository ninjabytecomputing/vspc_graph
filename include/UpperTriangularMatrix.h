#ifndef UPPER_TRIANGULAR_MATRIX_HAS_BEEN_INCLUDED
#define UPPER_TRIANGULAR_MATRIX_HAS_BEEN_INCLUDED

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

namespace vspc
{

template <typename T>
class UpperTriangularMatrix
{
public:
    using ValueType = T;

    explicit UpperTriangularMatrix(const size_t n = 0);

    UpperTriangularMatrix(const size_t n, const T& val);

    // Default copy, move, and destructor are fine

    /// @{
    /// @brief Two-dimensional index accessors.
    /// @details These are only provided for user-convenience. It is not
    /// recommended to use these for performant code.
    const T& operator()(const size_t i, const size_t j) const { return mData[index(i, j)]; }
          T& operator()(const size_t i, const size_t j)       { return mData[index(i, j)]; }
    /// @}

    /// @{
    /// @brief One-dimensional (flat) index accessors.
    /// @details These are the preferred way to access elements of this matrix
    /// class as long as the user has already precomputed the flat index (and
    /// is not recomputing the flat index on the fly).
    const T& operator()(const size_t idx) const { return mData[idx]; }
          T& operator()(const size_t idx)       { return mData[idx]; }
    /// @}

    /// @return Dimension of the square matrix.
    size_t dim() const { return mSize; }

    /// @return Number of elements in the upper triangular portion.
    /// @details This count also includes the diagonal elements.
    size_t numElements() const { return mData.size(); }

    /// @return Linear index into the upper-triangular matrix associated with
    /// the two-dimensional index.
    size_t index(size_t i, size_t j) const;

    /// @return Human-readable information about the matrix.
    std::string str() const;

private:
    size_t         mSize;
    std::vector<T> mData;
};

template <typename T>
std::ostream&
operator<<(std::ostream& os, const UpperTriangularMatrix<T>& obj)
{
    return os << obj.str();
}

// -----------------------------------------------------------------------------

template <typename T>
UpperTriangularMatrix<T>::UpperTriangularMatrix(const size_t n)
    : mSize(n), mData(n * (n + 1) / 2)
{
    // empty
}

template <typename T>
UpperTriangularMatrix<T>::UpperTriangularMatrix(const size_t n, const T& val)
    : mSize(n), mData(n * (n + 1) / 2, val)
{
    // empty
}

template <typename T>
size_t
UpperTriangularMatrix<T>::index(const size_t i, const size_t j) const
{
#ifdef _DEBUG
    assert(i <= j);
#endif
    return mData.size() - ((mSize - i) * (mSize + 1 - i) / 2 + i) + j;
}

template <typename T>
std::string
UpperTriangularMatrix<T>::str() const
{
    std::ostringstream os;
    os << "Matrix size: " << mSize << " x " << mSize << "\n";
    for (size_t i = 0; i < mSize; ++i) {
        for (size_t j = i; j < mSize; ++j) {
            os << mData[index(i, j)] << " ";
        }
        os << "\n";
    }
    os << "\n";
    return os.str();
}

} // namespace vspc

#endif // UPPER_TRIANGULAR_MATRIX_HAS_BEEN_INCLUDED
