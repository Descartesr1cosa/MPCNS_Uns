#include "induction/Faraday.h"

namespace mercury
{
using namespace uns;

std::vector<char> active_faces(const CellComplex &cx)
{
    std::vector<char> a(cx.n_faces(), 0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k) a[cx.c2f.idx[k]] = 1;
    return a;
}

void faraday_update(const CellComplex &cx, const std::vector<char> &active, const double *E, double *B, double dt)
{
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        if (!active[f]) continue;
        double circ = 0.0; // (d1 E)_f
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) circ += cx.f2e.sgn[k] * E[cx.f2e.idx[k]];
        B[f] -= dt * circ;
    }
}
} // namespace mercury
