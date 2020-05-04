#ifndef CSV_READER_HAS_BEEN_INCLUDED
#define CSV_READER_HAS_BEEN_INCLUDED

#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace vspc
{

class CSVRow
{
public:
    explicit CSVRow(const char delim = ',') : mDelim(delim) {}

    void readNextRow(std::istream& istr)
    {
        std::string line;
        std::getline(istr, line);

        std::istringstream lineStream(line);
        std::string cell;

        mData.clear();
        while(std::getline(lineStream, cell, mDelim)) {
            mData.push_back(cell);
        }

        // This checks for a trailing comma with no data after it.
        if (!lineStream && cell.empty()) {
            // If there was a trailing comma, then add an empty element.
            mData.push_back("");
        }
    }

    bool isMetadata() const
    {
        if (!mData.empty()) {
            const std::string& first = mData.front();
            if (!first.empty()) {
                return first.front() == '#';
            }
        }
        return false;
    }

    size_t size() const noexcept { return mData.size(); }

    char delimiter() const noexcept { return mDelim; }

    const std::string& operator[](size_t idx) const { return mData[idx]; }

private:
    const char               mDelim;
    std::vector<std::string> mData;
};

std::ostream& operator<<(std::ostream& ostr, const CSVRow& row)
{
    const size_t n = row.size();
    for (size_t i = 0; i < n-1; ++i) {
        ostr << row[i] << row.delimiter();
    }
    return ostr << row[n-1];
}

std::istream& operator>>(std::istream& istr, CSVRow& row)
{
    row.readNextRow(istr);
    return istr;
}

// -----------------------------------------------------------------------------

class CSVIterator
{
public:
    using iterator_category = std::input_iterator_tag;
    using value_type        = CSVRow;
    using difference_type   = size_t;
    using pointer           = CSVRow*;
    using reference         = CSVRow&;

    CSVIterator() : mStr(nullptr) {}

    CSVIterator(std::istream& str, const char delim = ',')
        : mStr(str.good() ? &str : nullptr), mRow(delim)
    {
        ++(*this);
    }

    CSVIterator& operator++()
    {
        if (mStr) {
            if (!((*mStr) >> mRow)) {
                mStr = nullptr;
            }
        }
        return *this;
    }

    CSVIterator operator++(int)
    {
        CSVIterator tmp(*this);
        ++(*this);
        return tmp;
    }

    operator bool() const { return mStr; }

    const CSVRow& operator*()  const { return mRow; }
    const CSVRow* operator->() const { return &mRow; }

private:
    std::istream* mStr;
    CSVRow        mRow;
};

} // namespace vspc

#endif // CSV_READER_HAS_BEEN_INCLUDED
