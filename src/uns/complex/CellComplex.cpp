#include "uns/complex/CellComplex.h"

#include <algorithm>
#include <stdexcept>

namespace uns
{
void CellComplex::build_up()
{
    n2e = e2n.transpose(n_nodes());
    e2f = f2e.transpose(n_edges());
    f2c = c2f.transpose(n_faces());
}

void CellComplex::edge_nodes(Id e, Id &start, Id &end) const
{
    start = end = -1;
    for (std::int64_t k = e2n.begin(e); k < e2n.end(e); ++k)
        (e2n.sgn[k] < 0 ? start : end) = e2n.idx[k];
}

void CellComplex::face_loop(Id f, std::vector<Id> &nodes) const
{
    nodes.clear();
    for (std::int64_t k = f2e.begin(f); k < f2e.end(f); ++k)
    {
        Id a, b;
        edge_nodes(f2e.idx[k], a, b);
        nodes.push_back(f2e.sgn[k] > 0 ? a : b); // node at which this edge is entered
    }
}

void CellComplex::face_loop_outward(Id, Id f, int sign, std::vector<Id> &nodes) const
{
    face_loop(f, nodes);
    // Reversing the traversal keeps the node set; the start node is irrelevant.
    if (sign < 0) std::reverse(nodes.begin(), nodes.end());
}

void CellComplex::cell_nodes(Id c, std::vector<Id> &nodes) const
{
    nodes.clear();
    std::vector<Id> loop;
    for (std::int64_t k = c2f.begin(c); k < c2f.end(c); ++k)
    {
        face_loop(c2f.idx[k], loop);
        nodes.insert(nodes.end(), loop.begin(), loop.end());
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
}
} // namespace uns

namespace uns
{
std::uint64_t topology_hash(const CellComplex &cx)
{
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&](const void *p, std::size_t n) {
        const unsigned char *b = static_cast<const unsigned char *>(p);
        for (std::size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
    };
    for (const Incidence *t : {&cx.e2n, &cx.f2e, &cx.c2f})
    {
        mix(t->ptr.data(), t->ptr.size() * sizeof(std::int64_t));
        mix(t->idx.data(), t->idx.size() * sizeof(Id));
        mix(t->sgn.data(), t->sgn.size());
    }
    return h;
}
} // namespace uns
