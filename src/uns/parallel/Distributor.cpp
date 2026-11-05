#include "uns/parallel/Distributor.h"

#include <algorithm>
#include <stdexcept>

namespace uns
{
Distributor::Distributor(const Mesh &global, const std::vector<int> &cell_rank, int nranks)
    : g_(global), nranks_(nranks)
{
    const CellComplex &cx = g_.cx;
    if (static_cast<Id>(cell_rank.size()) != cx.n_cells()) throw std::runtime_error("Distributor: partition size mismatch");
    for (int d = 0; d < kNumDims; ++d) owner_[d].assign(cx.count(static_cast<Dim>(d)), nranks);
    owner_[3] = cell_rank;
    // owner = min rank over the incident cells
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        const int p = cell_rank[c];
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            const Id f = cx.c2f.idx[kf];
            owner_[2][f] = std::min(owner_[2][f], p);
            for (auto ke = cx.f2e.begin(f); ke < cx.f2e.end(f); ++ke)
            {
                const Id e = cx.f2e.idx[ke];
                owner_[1][e] = std::min(owner_[1][e], p);
                for (auto kn = cx.e2n.begin(e); kn < cx.e2n.end(e); ++kn)
                    owner_[0][cx.e2n.idx[kn]] = std::min(owner_[0][cx.e2n.idx[kn]], p);
            }
        }
    }
    // node -> cells
    std::vector<std::int64_t> nptr(cx.n_nodes() + 1, 0);
    std::vector<Id> ncell;
    {
        std::vector<Id> cn;
        std::vector<std::pair<Id, Id>> pairs;
        for (Id c = 0; c < cx.n_cells(); ++c)
        {
            cx.cell_nodes(c, cn);
            for (Id v : cn) { pairs.emplace_back(v, c); ++nptr[v + 1]; }
        }
        for (Id v = 0; v < cx.n_nodes(); ++v) nptr[v + 1] += nptr[v];
        ncell.resize(pairs.size());
        std::vector<std::int64_t> fill(nptr.begin(), nptr.end() - 1);
        for (auto &pr : pairs) ncell[fill[pr.first]++] = pr.second;
    }

    local_.resize(nranks);
    n_owned_.resize(nranks);
    touch_edges_.resize(nranks);
    std::array<std::vector<int>, kNumDims> stamp;
    for (int d = 0; d < kNumDims; ++d) stamp[d].assign(cx.count(static_cast<Dim>(d)), -1);

    std::vector<std::vector<Id>> owned_cells(nranks);
    for (Id c = 0; c < cx.n_cells(); ++c) owned_cells[cell_rank[c]].push_back(c);

    for (int r = 0; r < nranks; ++r)
    {
        std::array<std::vector<GId>, kNumDims> found;
        auto mark = [&](int d, Id id) {
            if (stamp[d][id] != r) { stamp[d][id] = r; found[d].push_back(id); }
        };
        std::vector<Id> cn;
        for (Id c : owned_cells[r])
        {
            mark(3, c);
            cx.cell_nodes(c, cn);
            for (Id v : cn)
                for (auto k = nptr[v]; k < nptr[v + 1]; ++k) mark(3, ncell[k]);
        }
        for (GId c : found[3])
            for (auto kf = cx.c2f.begin(static_cast<Id>(c)); kf < cx.c2f.end(static_cast<Id>(c)); ++kf)
            {
                const Id f = cx.c2f.idx[kf];
                mark(2, f);
                for (auto ke = cx.f2e.begin(f); ke < cx.f2e.end(f); ++ke)
                {
                    const Id e = cx.f2e.idx[ke];
                    mark(1, e);
                    for (auto kn = cx.e2n.begin(e); kn < cx.e2n.end(e); ++kn) mark(0, cx.e2n.idx[kn]);
                }
            }
        for (int d = 0; d < kNumDims; ++d)
        {
            auto &v = found[d];
            const auto &own = owner_[d];
            std::sort(v.begin(), v.end(), [&](GId a, GId b) {
                const bool oa = own[a] == r, ob = own[b] == r;
                if (oa != ob) return oa;                       // owned first
                if (!oa && own[a] != own[b]) return own[a] < own[b];
                return a < b;
            });
            n_owned_[r][d] = static_cast<Id>(std::count_if(v.begin(), v.end(), [&](GId a) { return own[a] == r; }));
            local_[r][d] = std::move(v);
        }
        // edges touched by owned cells but owned by someone else
        std::vector<GId> touch;
        for (Id c : owned_cells[r])
            for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
            {
                const Id f = cx.c2f.idx[kf];
                for (auto ke = cx.f2e.begin(f); ke < cx.f2e.end(f); ++ke)
                    if (owner_[1][cx.f2e.idx[ke]] != r) touch.push_back(cx.f2e.idx[ke]);
            }
        std::sort(touch.begin(), touch.end());
        touch.erase(std::unique(touch.begin(), touch.end()), touch.end());
        std::sort(touch.begin(), touch.end(), [&](GId a, GId b) {
            return owner_[1][a] != owner_[1][b] ? owner_[1][a] < owner_[1][b] : a < b;
        });
        touch_edges_[r] = std::move(touch);
    }
}

LocalMesh Distributor::build(int r) const
{
    const CellComplex &gx = g_.cx;
    LocalMesh lm;
    lm.rank = r;
    lm.nranks = nranks_;
    CellComplex &lx = lm.mesh.cx;

    // global -> local maps
    std::array<std::vector<Id>, kNumDims> g2l;
    for (int d = 0; d < kNumDims; ++d)
    {
        g2l[d].assign(gx.count(static_cast<Dim>(d)), -1);
        for (std::size_t i = 0; i < local_[r][d].size(); ++i) g2l[d][local_[r][d][i]] = static_cast<Id>(i);
        EntityInfo &ei = lx.ent[d];
        ei.gid.assign(local_[r][d].begin(), local_[r][d].end());
        ei.n_owned = n_owned_[r][d];
        ei.owner.resize(ei.gid.size());
        for (std::size_t i = 0; i < ei.gid.size(); ++i) ei.owner[i] = owner_[d][ei.gid[i]];
        ei.alias_sign.assign(ei.gid.size(), 1); // canonical orientation derives from global node ids
    }
    auto copy_rows = [&](const Incidence &src, Incidence &dst, int row_dim, int col_dim) {
        for (GId gid : local_[r][row_dim])
        {
            const Id row = static_cast<Id>(gid);
            for (auto k = src.begin(row); k < src.end(row); ++k)
            {
                const Id l = g2l[col_dim][src.idx[k]];
                if (l < 0) throw std::runtime_error("Distributor: sub-complex is not closed");
                dst.push(l, src.sgn[k]);
            }
            dst.close_row();
        }
    };
    copy_rows(gx.e2n, lx.e2n, 1, 0);
    copy_rows(gx.f2e, lx.f2e, 2, 1);
    copy_rows(gx.c2f, lx.c2f, 3, 2);
    lx.build_up();
    lx.face_tag.resize(local_[r][2].size());
    for (std::size_t i = 0; i < local_[r][2].size(); ++i) lx.face_tag[i] = gx.face_tag[local_[r][2][i]];
    lx.cell_material.resize(local_[r][3].size());
    for (std::size_t i = 0; i < local_[r][3].size(); ++i) lx.cell_material[i] = gx.cell_material[local_[r][3][i]];
    lm.mesh.x.resize(local_[r][0].size());
    for (std::size_t i = 0; i < local_[r][0].size(); ++i) lm.mesh.x[i] = g_.x[local_[r][0][i]];
    lm.mesh.tag_names = g_.tag_names;

    // forward plans
    for (int d = 0; d < kNumDims; ++d)
    {
        const auto &own = owner_[d];
        // recv: my ghosts grouped by owner (already contiguous)
        {
            std::vector<Id> ids;
            int cur = -1;
            for (std::size_t i = n_owned_[r][d]; i < local_[r][d].size(); ++i)
            {
                const int o = own[local_[r][d][i]];
                if (o != cur)
                {
                    if (cur >= 0) lm.plan.recv[d].add_peer(cur, ids);
                    ids.clear();
                    cur = o;
                }
                ids.push_back(static_cast<Id>(i));
            }
            if (cur >= 0) lm.plan.recv[d].add_peer(cur, ids);
        }
        // send: for each peer q, ghosts of q that I own, in gid order (same order q uses to receive)
        for (int q = 0; q < nranks_; ++q)
        {
            if (q == r) continue;
            const auto &lq = local_[q][d];
            auto lo = std::lower_bound(lq.begin() + n_owned_[q][d], lq.end(), r, [&](GId a, int o) { return own[a] < o; });
            std::vector<Id> ids;
            for (auto it = lo; it != lq.end() && own[*it] == r; ++it) ids.push_back(g2l[d][*it]);
            if (!ids.empty()) lm.plan.send[d].add_peer(q, ids);
        }
    }
    // edge reduction plan: alias partial sums -> owner
    {
        const auto &own = owner_[1];
        const auto &te = touch_edges_[r];
        std::vector<Id> ids;
        int cur = -1;
        for (GId e : te)
        {
            if (own[e] != cur)
            {
                if (cur >= 0) lm.plan.reduce_send.add_peer(cur, ids);
                ids.clear();
                cur = own[e];
            }
            ids.push_back(g2l[1][e]);
        }
        if (cur >= 0) lm.plan.reduce_send.add_peer(cur, ids);
        for (int q = 0; q < nranks_; ++q)
        {
            if (q == r) continue;
            const auto &tq = touch_edges_[q];
            auto lo = std::lower_bound(tq.begin(), tq.end(), r, [&](GId a, int o) { return own[a] < o; });
            std::vector<Id> rid;
            for (auto it = lo; it != tq.end() && own[*it] == r; ++it) rid.push_back(g2l[1][*it]);
            if (!rid.empty()) lm.plan.reduce_recv.add_peer(q, rid);
        }
    }
    return lm;
}
} // namespace uns
