#include "uns/complex/Incidence.h"

#include <algorithm>
#include <map>

namespace uns
{
Incidence Incidence::transpose(Id ncols) const
{
    Incidence t;
    t.ptr.assign(static_cast<std::size_t>(ncols) + 1, 0);
    for (Id c : idx) ++t.ptr[static_cast<std::size_t>(c) + 1];
    for (Id c = 0; c < ncols; ++c) t.ptr[c + 1] += t.ptr[c];
    t.idx.resize(idx.size());
    t.sgn.resize(sgn.size());
    std::vector<std::int64_t> fill(t.ptr.begin(), t.ptr.end() - 1);
    // Rows are visited in ascending order, so every new row is sorted by column.
    for (Id r = 0; r < rows(); ++r)
        for (std::int64_t k = ptr[r]; k < ptr[r + 1]; ++k)
        {
            const std::int64_t p = fill[idx[k]]++;
            t.idx[p] = r;
            t.sgn[p] = sgn[k];
        }
    return t;
}

long count_nonzero_products(const Incidence &upper, const Incidence &lower)
{
    long bad = 0;
    std::map<Id, long> acc;
    for (Id r = 0; r < upper.rows(); ++r)
    {
        acc.clear();
        for (std::int64_t k = upper.begin(r); k < upper.end(r); ++k)
        {
            const Id m = upper.idx[k];
            for (std::int64_t q = lower.begin(m); q < lower.end(m); ++q)
                acc[lower.idx[q]] += static_cast<long>(upper.sgn[k]) * lower.sgn[q];
        }
        for (const auto &kv : acc)
            if (kv.second != 0) { ++bad; break; }
    }
    return bad;
}
} // namespace uns
