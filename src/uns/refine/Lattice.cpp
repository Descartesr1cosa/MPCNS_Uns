#include "refine/Lattice.h"

#include <algorithm>
#include <set>
#include <stdexcept>

namespace uns::refine_detail
{
namespace
{
std::uint64_t pair_key(Id a, Id b)
{
    const Id lo = std::min(a, b), hi = std::max(a, b);
    return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
}
} // namespace

Lattice::Lattice(const Mesh &m, int max_level, bool snap) : m_(m), Lmax_(max_level), N_(1 << max_level), snap_(snap)
{
    const CellComplex &cx = m.cx;
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        Id a, b;
        cx.edge_nodes(e, a, b);
        edge_of_[pair_key(a, b)] = e;
    }
    std::map<std::array<Id, 4>, Id> face_of;
    std::vector<Id> loop;
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        cx.face_loop(f, loop);
        if (loop.size() != 4) continue;
        std::array<Id, 4> k = {loop[0], loop[1], loop[2], loop[3]};
        std::sort(k.begin(), k.end());
        face_of[k] = f;
    }
    fr_.resize(cx.n_cells());
    std::vector<Id> nodes;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        cx.cell_nodes(c, nodes);
        if (nodes.size() != 8 || cx.c2f.count(c) != 6)
            throw std::runtime_error("refine: coarse cell " + std::to_string(cx.ent[3].gid[c]) +
                                     " is not a hexahedron (only hexahedral coarse cells can be refined)");
        std::map<Id, std::set<Id>> nb;
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            cx.face_loop(cx.c2f.idx[kf], loop);
            if (loop.size() != 4) throw std::runtime_error("refine: non-quadrilateral face on a coarse hex");
            for (int i = 0; i < 4; ++i)
            {
                nb[loop[i]].insert(loop[(i + 1) % 4]);
                nb[loop[(i + 1) % 4]].insert(loop[i]);
            }
        }
        for (Id n : nodes)
            if (nb[n].size() != 3) throw std::runtime_error("refine: coarse cell is not a topological hexahedron");
        CellFrame &f = fr_[c];
        const Id n0 = nodes[0];
        auto it = nb[n0].begin();
        Id n1 = *it++, n2 = *it++, n3 = *it;
        if (dot(cross(m.x[n1] - m.x[n0], m.x[n2] - m.x[n0]), m.x[n3] - m.x[n0]) < 0) std::swap(n2, n3);
        auto common = [&](Id a, Id b) {
            Id r = -1;
            for (Id x : nb[a])
                if (x != n0 && nb[b].count(x)) r = x;
            if (r < 0) throw std::runtime_error("refine: cannot identify hex corner");
            return r;
        };
        f.node[0] = n0; f.node[1] = n1; f.node[2] = n2; f.node[4] = n3;
        f.node[3] = common(n1, n2); f.node[5] = common(n1, n3); f.node[6] = common(n2, n3);
        f.node[7] = -1;
        for (Id n : nodes)
            if (std::find(f.node, f.node + 7, n) == f.node + 7) f.node[7] = n;
        for (int ax = 0; ax < 3; ++ax)
            for (int pos = 0; pos < 2; ++pos)
            {
                const int u = (ax + 1) % 3, v = (ax + 2) % 3;
                Id corner[2][2];
                std::array<Id, 4> k;
                for (int iu = 0; iu < 2; ++iu)
                    for (int iv = 0; iv < 2; ++iv)
                        corner[iu][iv] = f.node[(pos << ax) | (iu << u) | (iv << v)];
                k = {corner[0][0], corner[0][1], corner[1][0], corner[1][1]};
                std::sort(k.begin(), k.end());
                auto fit = face_of.find(k);
                if (fit == face_of.end()) throw std::runtime_error("refine: hex side face not found");
                SideFrame &sf = f.side[2 * ax + pos];
                sf.face = fit->second;
                cx.face_loop(sf.face, loop);
                auto where = [&](Id n, int &iu, int &iv) {
                    for (iu = 0; iu < 2; ++iu)
                        for (iv = 0; iv < 2; ++iv)
                            if (corner[iu][iv] == n) return;
                    throw std::runtime_error("refine: face loop node not on hex side");
                };
                int u0, v0, u1, v1, u3, v3;
                where(loop[0], u0, v0); where(loop[1], u1, v1); where(loop[3], u3, v3);
                sf.ou = u0; sf.ov = v0;
                sf.sa = (u1 != u0) ? 0 : 1;
                sf.ta = (u3 != u0) ? 0 : 1;
                if (sf.sa == sf.ta) throw std::runtime_error("refine: degenerate face frame");
            }
    }
}

void Lattice::side_st(Id cell, int side, int cu, int cv, int &s, int &t) const
{
    const SideFrame &sf = fr_[cell].side[side];
    const int c[2] = {cu, cv};
    const int o[2] = {sf.ou * N_, sf.ov * N_};
    s = std::abs(c[sf.sa] - o[sf.sa]);
    t = std::abs(c[sf.ta] - o[sf.ta]);
}

NodeKey Lattice::key(Id cell, int a, int b, int c) const
{
    const int p[3] = {a, b, c};
    int nb = 0;
    for (int k = 0; k < 3; ++k) nb += (p[k] == 0 || p[k] == N_);
    const CellFrame &f = fr_[cell];
    NodeKey k;
    if (nb == 3)
    {
        const int q = (p[0] == N_ ? 1 : 0) | (p[1] == N_ ? 2 : 0) | (p[2] == N_ ? 4 : 0);
        k.a = (0ull << 56) | static_cast<std::uint32_t>(f.node[q]);
    }
    else if (nb == 2)
    {
        int ax = 0;
        for (int i = 0; i < 3; ++i)
            if (p[i] != 0 && p[i] != N_) ax = i;
        const int q0 = (p[0] == N_ ? 1 : 0) | (p[1] == N_ ? 2 : 0) | (p[2] == N_ ? 4 : 0);
        const Id n0 = f.node[q0 & ~(1 << ax)], n1 = f.node[q0 | (1 << ax)];
        const Id e = edge_of_.at(pair_key(n0, n1));
        Id start, end;
        m_.cx.edge_nodes(e, start, end);
        const std::uint64_t pos = (n0 == start) ? p[ax] : N_ - p[ax];
        k.a = (1ull << 56) | static_cast<std::uint32_t>(e);
        k.b = pos;
    }
    else if (nb == 1)
    {
        int ax = 0;
        for (int i = 0; i < 3; ++i)
            if (p[i] == 0 || p[i] == N_) ax = i;
        const int side = 2 * ax + (p[ax] == N_ ? 1 : 0);
        int s, t;
        side_st(cell, side, p[(ax + 1) % 3], p[(ax + 2) % 3], s, t);
        k.a = (2ull << 56) | static_cast<std::uint32_t>(f.side[side].face);
        k.b = static_cast<std::uint64_t>(s) | (static_cast<std::uint64_t>(t) << 21);
    }
    else
    {
        k.a = (3ull << 56) | static_cast<std::uint32_t>(cell);
        k.b = static_cast<std::uint64_t>(a) | (static_cast<std::uint64_t>(b) << 21) | (static_cast<std::uint64_t>(c) << 42);
    }
    return k;
}

Vec3 Lattice::position(Id cell, int a, int b, int c) const
{
    const CellFrame &f = fr_[cell];
    const double x[3] = {double(a) / N_, double(b) / N_, double(c) / N_};
    Vec3 X;
    double rad = 0;
    for (int q = 0; q < 8; ++q)
    {
        const double w = ((q & 1) ? x[0] : 1 - x[0]) * ((q & 2) ? x[1] : 1 - x[1]) * ((q & 4) ? x[2] : 1 - x[2]);
        X += w * m_.x[f.node[q]];
        rad += w * norm(m_.x[f.node[q]]);
    }
    if (snap_)
    {
        const double n = norm(X);
        if (n > 0) X *= rad / n;
    }
    return X;
}
} // namespace uns::refine_detail
