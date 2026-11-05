#pragma once
// Edge assembly shared by every EMF term.
//
//   begin(E)                       zero the edge array
//   add_cell_edge_term(c, e, v)    a cell-based contribution (already weighted by w_{e,c}
//                                  and projected on the edge: line-integral units)
//   finish(E)                      owner reduction + broadcast  (Exchanger::assemble_edges)
//
// Ideal, Hall, resistive and artificial terms all call add_cell_edge_term with their own
// values; the weights w_{e,c} come from uns::Geometry (dual-face shares), so the same code
// serves any edge valence, hanging nodes and MPI interfaces.
#include <algorithm>

#include "uns/comm/Exchanger.h"
#include "uns/geometry/Geometry.h"

namespace mercury
{
class EdgeEMF
{
public:
    EdgeEMF(const uns::LocalMesh &lm, const uns::Geometry &g, uns::comm::Exchanger &ex) : lm_(lm), g_(g), ex_(ex) {}
    void begin(double *E) const;
    // weight w_{e,c} of owned cell c at edge e (looked up in the geometry CSR)
    double weight(uns::Id e, uns::Id c) const;
    void add(double *E, uns::Id e, double v) const { E[e] += v; }
    void finish(double *E) { ex_.assemble_edges(E); }
    // visit (c, e, w) for every owned cell and every edge of it
    template <class F> void for_each_owned_cell_edge(F &&f) const;
private:
    const uns::LocalMesh &lm_;
    const uns::Geometry &g_;
    uns::comm::Exchanger &ex_;
};

template <class F> void EdgeEMF::for_each_owned_cell_edge(F &&f) const
{
    const uns::CellComplex &cx = lm_.mesh.cx;
    std::vector<uns::Id> es;
    for (uns::Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        es.clear();
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
            for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke) es.push_back(cx.f2e.idx[ke]);
        std::sort(es.begin(), es.end());
        es.erase(std::unique(es.begin(), es.end()), es.end());
        for (uns::Id e : es) f(c, e, weight(e, c));
    }
}
} // namespace mercury
