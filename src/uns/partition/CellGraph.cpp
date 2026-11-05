#include "uns/partition/CellGraph.h"

#include <algorithm>

namespace uns
{
CellGraph build_cell_graph(const Mesh &m, const Geometry &g, const CostModel &cost)
{
    const CellComplex &cx = m.cx;
    const Id nc = cx.n_cells();
    CellGraph out;
    out.vw.resize(nc);
    out.centre = g.cell_center;
    for (Id c = 0; c < nc; ++c)
    {
        const double base = cx.cell_material[c] == kSolid ? cost.solid : cost.fluid;
        out.vw[c] = base * (1.0 + cost.per_extra_face * std::max(0, cx.c2f.count(c) - 6));
    }
    out.xadj.assign(nc + 1, 0);
    // two passes over interior faces
    std::vector<int> deg(nc, 0);
    for (Id f = 0; f < cx.n_faces(); ++f)
        if (cx.f2c.count(f) == 2)
        {
            ++deg[cx.f2c.idx[cx.f2c.begin(f)]];
            ++deg[cx.f2c.idx[cx.f2c.begin(f) + 1]];
        }
    for (Id c = 0; c < nc; ++c) out.xadj[c + 1] = out.xadj[c] + deg[c];
    out.adj.resize(out.xadj[nc]);
    out.adjw.resize(out.xadj[nc]);
    std::vector<std::int64_t> fill(out.xadj.begin(), out.xadj.end() - 1);
    for (Id f = 0; f < cx.n_faces(); ++f)
        if (cx.f2c.count(f) == 2)
        {
            const Id a = cx.f2c.idx[cx.f2c.begin(f)], b = cx.f2c.idx[cx.f2c.begin(f) + 1];
            const bool solid = cx.cell_material[a] == kSolid || cx.cell_material[b] == kSolid;
            const double w = (solid ? cost.comm_solid : cost.comm_fluid);
            out.adj[fill[a]] = b; out.adjw[fill[a]++] = w;
            out.adj[fill[b]] = a; out.adjw[fill[b]++] = w;
        }
    return out;
}
} // namespace uns
