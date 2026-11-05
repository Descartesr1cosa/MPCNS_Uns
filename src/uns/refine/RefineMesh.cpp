// Octree refinement of hexahedral coarse cells into a non-conforming polyhedral complex.
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <unordered_map>

#include "refine/Lattice.h"
#include "uns/complex/ComplexBuilder.h"
#include "uns/refine/Refine.h"

namespace uns
{
using namespace refine_detail;

namespace
{
struct Leaf
{
    Id cell;
    int level;
    int p[3]; // anchor in finest-lattice units
};

struct Entry // one side-square of one leaf
{
    std::uint64_t pk_a, pk_b; // plane key
    int grp;                  // which side of the plane
    int level, u, v;          // dyadic square in plane coordinates
    int leaf, side;
    bool matched = false, covered = false;
};

struct FaceRec
{
    int owner_leaf, owner_side, other_leaf; // other_leaf = -1: domain boundary
    Id coarse_face;                         // -1: interior of a coarse cell
};

inline std::uint64_t sq_key(int level, int u, int v)
{
    return (static_cast<std::uint64_t>(level) << 56) | (static_cast<std::uint64_t>(u) << 28) | static_cast<std::uint64_t>(v);
}

struct Analysis
{
    std::vector<FaceRec> faces;
    std::vector<char> split; // leaves that must be split to satisfy the face balance
    bool any_split = false;
};

class Builder
{
public:
    Builder(const Mesh &coarse, const RefineSpec &spec, const std::vector<int> &cell_levels, const RefineOptions &opt)
        : cx_(coarse.cx), spec_(spec), levels_(cell_levels), opt_(opt), lat_(coarse, spec.max_level, opt.spherical_snap)
    {
        N_ = lat_.N();
        Lmax_ = lat_.max_level();
    }

    RefineResult run(const Mesh &coarse);

private:
    int wanted_level(const Leaf &l) const;
    void generate(Leaf l, std::vector<Leaf> &out) const;
    void children(const Leaf &l, std::vector<Leaf> &out) const;
    Analysis analyse(const std::vector<Leaf> &leaves) const;
    void corners_of(const Leaf &l, int side, int pt[4][3]) const;
    void segment(Id cell, const int P[3], const int Q[3], std::vector<Id> &out,
                 const std::unordered_map<NodeKey, Id, NodeKeyHash> &ids) const;

    const CellComplex &cx_;
    const RefineSpec &spec_;
    const std::vector<int> &levels_;
    RefineOptions opt_;
    Lattice lat_;
    int N_, Lmax_;
};

int Builder::wanted_level(const Leaf &l) const
{
    int w = levels_[l.cell];
    const int s = N_ >> l.level;
    if (!spec_.regions.empty())
    {
        int reg = spec_.region_level(lat_.position(l.cell, l.p[0] + s / 2, l.p[1] + s / 2, l.p[2] + s / 2));
        for (int q = 0; q < 8 && reg < spec_.max_level; ++q)
            reg = std::max(reg, spec_.region_level(lat_.position(l.cell, l.p[0] + ((q & 1) ? s : 0), l.p[1] + ((q & 2) ? s : 0),
                                                                  l.p[2] + ((q & 4) ? s : 0))));
        w = std::max(w, reg);
    }
    if (cx_.cell_material[l.cell] == kSolid) w = std::min(w, spec_.solid_max_level);
    return std::min(w, Lmax_);
}

void Builder::children(const Leaf &l, std::vector<Leaf> &out) const
{
    const int h = (N_ >> l.level) / 2;
    for (int q = 0; q < 8; ++q)
        out.push_back({l.cell, l.level + 1, {l.p[0] + ((q & 1) ? h : 0), l.p[1] + ((q & 2) ? h : 0), l.p[2] + ((q & 4) ? h : 0)}});
}

void Builder::generate(Leaf l, std::vector<Leaf> &out) const
{
    if (l.level < Lmax_ && l.level < wanted_level(l))
    {
        std::vector<Leaf> ch;
        children(l, ch);
        for (const Leaf &c : ch) generate(c, out);
    }
    else
        out.push_back(l);
}

void Builder::corners_of(const Leaf &l, int side, int pt[4][3]) const
{
    const int ax = side / 2, hi = side & 1, s = N_ >> l.level;
    const int u = (ax + 1) % 3, v = (ax + 2) % 3;
    static const int up[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};   // counter-clockwise seen from +axis
    static const int dn[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};   // outward for the low side
    for (int k = 0; k < 4; ++k)
    {
        const int *o = hi ? up[k] : dn[k];
        pt[k][ax] = l.p[ax] + (hi ? s : 0);
        pt[k][u] = l.p[u] + o[0] * s;
        pt[k][v] = l.p[v] + o[1] * s;
    }
}

Analysis Builder::analyse(const std::vector<Leaf> &leaves) const
{
    std::vector<Entry> ent;
    ent.reserve(leaves.size() * 6);
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
        const Leaf &l = leaves[i];
        const int s = N_ >> l.level;
        for (int sd = 0; sd < 6; ++sd)
        {
            const int ax = sd / 2, hi = sd & 1;
            const int pos = l.p[ax] + (hi ? s : 0);
            const int u = (ax + 1) % 3, v = (ax + 2) % 3;
            Entry e;
            e.level = l.level; e.leaf = static_cast<int>(i); e.side = sd;
            if ((hi && pos == N_) || (!hi && pos == 0))
            {
                const Id F = lat_.frame(l.cell).side[sd].face;
                int s0, t0, s1, t1;
                lat_.side_st(l.cell, sd, l.p[u], l.p[v], s0, t0);
                lat_.side_st(l.cell, sd, l.p[u] + s, l.p[v] + s, s1, t1);
                e.pk_a = (0ull << 60) | static_cast<std::uint32_t>(F);
                e.pk_b = 0;
                e.u = std::min(s0, s1);
                e.v = std::min(t0, t1);
                e.grp = (l.cell != cx_.f2c.idx[cx_.f2c.begin(F)]) ? 1 : 0;
            }
            else
            {
                e.pk_a = (1ull << 60) | static_cast<std::uint32_t>(l.cell);
                e.pk_b = (static_cast<std::uint64_t>(ax) << 40) | static_cast<std::uint64_t>(pos);
                e.u = l.p[u];
                e.v = l.p[v];
                e.grp = hi;
            }
            ent.push_back(e);
        }
    }
    std::vector<int> order(ent.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    std::sort(order.begin(), order.end(), [&](int x, int y) {
        const Entry &a = ent[x], &b = ent[y];
        if (a.pk_a != b.pk_a) return a.pk_a < b.pk_a;
        if (a.pk_b != b.pk_b) return a.pk_b < b.pk_b;
        if (a.grp != b.grp) return a.grp < b.grp;
        return x < y;
    });

    Analysis an;
    an.split.assign(leaves.size(), 0);
    auto record = [&](const Entry &owner, int other_leaf, const Entry &big_or_self) {
        FaceRec f;
        f.owner_leaf = owner.leaf; f.owner_side = owner.side; f.other_leaf = other_leaf;
        f.coarse_face = (owner.pk_a >> 60) == 0 ? static_cast<Id>(owner.pk_a & 0xffffffffu) : -1;
        an.faces.push_back(f);
        if (other_leaf >= 0 && big_or_self.level + 1 < owner.level)
        {
            an.split[other_leaf] = 1;
            an.any_split = true;
        }
    };
    std::unordered_map<std::uint64_t, int> hP, hM;
    for (std::size_t g0 = 0; g0 < order.size();)
    {
        std::size_t g1 = g0;
        while (g1 < order.size() && ent[order[g1]].pk_a == ent[order[g0]].pk_a && ent[order[g1]].pk_b == ent[order[g0]].pk_b) ++g1;
        hP.clear(); hM.clear();
        for (std::size_t k = g0; k < g1; ++k)
        {
            const Entry &e = ent[order[k]];
            (e.grp == 0 ? hP : hM)[sq_key(e.level, e.u, e.v)] = order[k];
        }
        auto find_cover = [&](const std::unordered_map<std::uint64_t, int> &h, const Entry &e, bool allow_equal) -> int {
            for (int lev = allow_equal ? e.level : e.level - 1; lev >= 0; --lev)
            {
                const int sz = N_ >> lev;
                const int u = e.u & ~(sz - 1), v = e.v & ~(sz - 1);
                auto it = h.find(sq_key(lev, u, v));
                if (it != h.end()) return it->second;
            }
            return -1;
        };
        for (std::size_t k = g0; k < g1; ++k)
        {
            Entry &e = ent[order[k]];
            if (e.grp == 0)
            {
                const int hit = find_cover(hM, e, true);
                if (hit >= 0)
                {
                    Entry &o = ent[hit];
                    e.matched = true;
                    o.matched |= (o.level == e.level);
                    o.covered |= (o.level < e.level);
                    record(e, o.leaf, o);
                }
            }
            else
            {
                const int hit = find_cover(hP, e, false);
                if (hit >= 0)
                {
                    Entry &o = ent[hit];
                    e.matched = true;
                    o.covered = true;
                    record(e, o.leaf, o);
                }
            }
        }
        for (std::size_t k = g0; k < g1; ++k)
        {
            Entry &e = ent[order[k]];
            if (!e.matched && !e.covered) record(e, -1, e);
        }
        g0 = g1;
    }
    return an;
}

void Builder::segment(Id cell, const int P[3], const int Q[3], std::vector<Id> &out,
                      const std::unordered_map<NodeKey, Id, NodeKeyHash> &ids) const
{
    int len = 0;
    for (int k = 0; k < 3; ++k) len = std::max(len, std::abs(Q[k] - P[k]));
    if (len < 2) return;
    int M[3];
    for (int k = 0; k < 3; ++k) M[k] = (P[k] + Q[k]) / 2;
    auto it = ids.find(lat_.key(cell, M[0], M[1], M[2]));
    if (it == ids.end()) return;
    segment(cell, P, M, out, ids);
    out.push_back(it->second);
    segment(cell, M, Q, out, ids);
}

RefineResult Builder::run(const Mesh &coarse)
{
    if (static_cast<Id>(levels_.size()) != cx_.n_cells()) throw std::runtime_error("refine: cell level array size mismatch");
    std::vector<Leaf> leaves;
    for (Id c = 0; c < cx_.n_cells(); ++c) generate({c, 0, {0, 0, 0}}, leaves);

    RefineResult rr;
    Analysis an;
    for (int it = 0;; ++it)
    {
        an = analyse(leaves);
        if (!an.any_split || it >= opt_.max_balance_iter) break;
        std::vector<Leaf> next;
        next.reserve(leaves.size() + 8);
        for (std::size_t i = 0; i < leaves.size(); ++i)
            if (an.split[i]) children(leaves[i], next); else next.push_back(leaves[i]);
        leaves.swap(next);
        rr.balance_iterations = it + 1;
        if (leaves.size() > 400000000u) throw std::runtime_error("refine: more than 4e8 leaves, aborting");
    }
    if (an.any_split) throw std::runtime_error("refine: 2:1 balance did not converge");

    // nodes: every leaf corner, in leaf order
    std::unordered_map<NodeKey, Id, NodeKeyHash> ids;
    std::vector<Vec3> xs;
    ids.reserve(leaves.size() * 2);
    for (const Leaf &l : leaves)
    {
        const int s = N_ >> l.level;
        for (int q = 0; q < 8; ++q)
        {
            const int a = l.p[0] + ((q & 1) ? s : 0), b = l.p[1] + ((q & 2) ? s : 0), c = l.p[2] + ((q & 4) ? s : 0);
            const NodeKey k = lat_.key(l.cell, a, b, c);
            if (ids.emplace(k, static_cast<Id>(xs.size())).second) xs.push_back(lat_.position(l.cell, a, b, c));
        }
    }
    // face loops (with hanging nodes) and per-leaf face lists
    std::vector<std::vector<Id>> loops(an.faces.size());
    std::vector<int> cnt(leaves.size() + 1, 0);
    for (std::size_t f = 0; f < an.faces.size(); ++f)
    {
        const FaceRec &fr = an.faces[f];
        const Leaf &l = leaves[fr.owner_leaf];
        int pt[4][3];
        corners_of(l, fr.owner_side, pt);
        for (int k = 0; k < 4; ++k)
        {
            loops[f].push_back(ids.at(lat_.key(l.cell, pt[k][0], pt[k][1], pt[k][2])));
            segment(l.cell, pt[k], pt[(k + 1) % 4], loops[f], ids);
        }
        ++cnt[fr.owner_leaf + 1];
        if (fr.other_leaf >= 0) ++cnt[fr.other_leaf + 1];
    }
    for (std::size_t i = 0; i < leaves.size(); ++i) cnt[i + 1] += cnt[i];
    std::vector<std::pair<int, bool>> lf(cnt.back()); // (face rec, reversed)
    std::vector<int> fill(cnt.begin(), cnt.end() - 1);
    for (std::size_t f = 0; f < an.faces.size(); ++f)
    {
        lf[fill[an.faces[f].owner_leaf]++] = {static_cast<int>(f), false};
        if (an.faces[f].other_leaf >= 0) lf[fill[an.faces[f].other_leaf]++] = {static_cast<int>(f), true};
    }

    ComplexBuilder b(static_cast<Id>(xs.size()));
    std::vector<std::vector<Id>> cell_loops;
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
        cell_loops.clear();
        for (int k = cnt[i]; k < cnt[i + 1]; ++k)
        {
            cell_loops.push_back(loops[lf[k].first]);
            if (lf[k].second) std::reverse(cell_loops.back().begin(), cell_loops.back().end());
        }
        b.add_cell(cell_loops, cx_.cell_material[leaves[i].cell]);
    }
    for (std::size_t f = 0; f < an.faces.size(); ++f)
        if (an.faces[f].coarse_face >= 0 && cx_.face_tag[an.faces[f].coarse_face] != 0)
            b.tag_face(loops[f], cx_.face_tag[an.faces[f].coarse_face]);

    rr.fine.cx = b.finish();
    rr.fine.x = std::move(xs);
    rr.fine.tag_names = coarse.tag_names;
    rr.cell_level.resize(leaves.size());
    rr.parent_cell.resize(leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
        rr.cell_level[i] = leaves[i].level;
        rr.parent_cell[i] = leaves[i].cell;
    }
    rr.face_coarse_face.assign(rr.fine.cx.n_faces(), -1);
    rr.face_owner_cell.assign(rr.fine.cx.n_faces(), -1);
    // loops were registered in creation order of faces inside ComplexBuilder; look them up by node set
    for (std::size_t f = 0; f < an.faces.size(); ++f)
    {
        const Id fid = b.face_id(loops[f]);
        rr.face_coarse_face[fid] = an.faces[f].coarse_face;
        rr.face_owner_cell[fid] = an.faces[f].owner_leaf;
    }
    return rr;
}
} // namespace

RefineResult build_refined_mesh(const Mesh &coarse, const RefineSpec &spec, const std::vector<int> &cell_levels,
                                const RefineOptions &opt)
{
    Builder b(coarse, spec, cell_levels, opt);
    return b.run(coarse);
}
} // namespace uns
