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

    size_t dim() const { return mSize; }

    size_t index(size_t i, size_t j) const;

    const T& operator()(const size_t i, const size_t j) const { return mData[index(i, j)]; }
          T& operator()(const size_t i, const size_t j)       { return mData[index(i, j)]; }

    const T& operator()(const size_t idx) const { return mData[idx]; }
          T& operator()(const size_t idx)       { return mData[idx]; }

    std::string str() const;

private:
    size_t         mSize;
    std::vector<T> mData;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const UpperTriangularMatrix<T>& obj)
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
    : UpperTriangularMatrix(n)
{
    std::fill(mData.begin(), mData.end(), val);
}

template <typename T>
size_t
UpperTriangularMatrix<T>::index(const size_t i, const size_t j) const
{
    assert(i <= j);
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
