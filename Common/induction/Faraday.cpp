#include "induction/Faraday.h"

namespace mhd
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

// B -= d1(dt_e E): the cochain dt_e E_e is what enters d1, so d2 B is preserved for ANY per-edge weights (pseudo-time stepping)
void faraday_update_edge_dt(const CellComplex &cx, const std::vector<char> &active, const double *E, const double *dt_edge, double *B)
{
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        if (!active[f]) continue;
        double circ = 0.0;
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) circ += cx.f2e.sgn[k] * dt_edge[cx.f2e.idx[k]] * E[cx.f2e.idx[k]];
        B[f] -= circ;
    }
}
} // namespace mhd
