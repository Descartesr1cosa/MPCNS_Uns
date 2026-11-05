#include "diagnostics/DivB.h"

#include <algorithm>
#include <cmath>

namespace mercury
{
DivBReport div_b(const uns::CellComplex &cx, const double *B, const double *extra)
{
    DivBReport r;
    for (uns::Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        double s = 0, a = 0;
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            s += cx.c2f.sgn[k] * B[cx.c2f.idx[k]];
            a += std::abs(B[cx.c2f.idx[k]]) + (extra ? std::abs(extra[cx.c2f.idx[k]]) : 0.0);
        }
        const double rel = a > 0 ? std::abs(s) / a : 0.0;
        if (std::abs(s) > r.max_abs) { r.max_abs = std::abs(s); r.worst_cell_gid = cx.ent[3].gid[c]; }
        r.max_rel = std::max(r.max_rel, rel);
    }
    return r;
}
} // namespace mercury
