#pragma once
// Signed CSR incidence table: row -> (column, sign).
#include <cstdint>
#include <vector>

#include "uns/core/Types.h"

namespace uns
{
struct Incidence
{
    std::vector<std::int64_t> ptr{0};
    std::vector<Id> idx;
    std::vector<std::int8_t> sgn;

    Id rows() const { return static_cast<Id>(ptr.size()) - 1; }
    std::int64_t nnz() const { return static_cast<std::int64_t>(idx.size()); }
    std::int64_t begin(Id r) const { return ptr[r]; }
    std::int64_t end(Id r) const { return ptr[r + 1]; }
    int count(Id r) const { return static_cast<int>(ptr[r + 1] - ptr[r]); }

    void clear() { ptr.assign(1, 0); idx.clear(); sgn.clear(); }
    void push(Id col, int sign) { idx.push_back(col); sgn.push_back(static_cast<std::int8_t>(sign)); }
    void close_row() { ptr.push_back(static_cast<std::int64_t>(idx.size())); }

    // Transposed table; entries of each new row are sorted by ascending column.
    Incidence transpose(Id ncols) const;
};

// Exact integer check of the composition  upper * lower == 0
// (upper: rows -> middle ids, lower: middle ids -> columns).
// Returns the number of rows whose composition is not identically zero.
long count_nonzero_products(const Incidence &upper, const Incidence &lower);
} // namespace uns
