#include "multigrid/Multigrid.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <cstdio>
#include <numeric>
#include <tuple>

#include "fluid/FluidRHS.h"
#include "physics/Euler.h"
#include "physics/Gyration.h"

namespace mercury
{
using namespace uns;

double positivity_safeguard(const double *U, const double *dU, double gamma, double keep, double grow)
{
    (void)gamma;
    auto internal = [](const double *u) { return u[4] - (u[0] > 0 ? 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0] : 0.0); };
    const double e_old = internal(U);
    double a = 1.0;
    for (int it = 0; it < 8; ++it, a *= 0.5)
    {
        double u[5];
        for (int m = 0; m < 5; ++m) u[m] = U[m] + a * dU[m];
        if (!(u[0] >= keep * U[0])) continue;
        if (e_old > 0 && !(internal(u) >= keep * e_old)) continue;
        if (!(internal(u) > 0)) continue;
        if (grow > 0 && (u[0] > grow * U[0] || (e_old > 0 && internal(u) > grow * e_old))) continue;
        return a;
    }
    return 0.0;
}

Multigrid::Multigrid(MPI_Comm comm, const LocalMesh &lm, const Geometry &geo, const BoundaryInfo &bc, comm::Exchanger &ex, const Config &cfg, const Prim fs[2])
    : comm_(comm), lm_(lm), geo_(geo), bc_(bc), ex_(ex), cfg_(cfg)
{
    fs_[0] = fs[0]; fs_[1] = fs[1];
    MPI_Comm_size(comm, &nrank_);
    const CellComplex &cx = lm.mesh.cx;
    const Id nc = cx.n_cells(), no = cx.ent[3].n_owned;
    const int nlev = std::max(1, cfg.mg_levels);
    lv_.resize(nlev);
    std::vector<Id> own_id(nc, -1);
    Id n_own = 0;
    for (Id c = 0; c < no; ++c)
        if (cx.cell_material[c] != kSolid) own_id[c] = n_own++;
    for (int l = 0; l < nlev; ++l)
    {
        build_level(l, own_id, n_own);
        if (l + 1 < nlev)
        {
            std::vector<Id> parent;
            const Id n_new = match(l, parent);
            lv_[l].parent = parent;
            for (Id c = 0; c < no; ++c)
                if (own_id[c] >= 0) own_id[c] = parent[own_id[c]];
            n_own = n_new;
        }
    }
    for (MgLevel &L : lv_)
    {
        const std::size_t nt = L.n_tot(), n = L.n_own;
        for (int k = 0; k < 2; ++k)
        {
            L.U[k].assign(5 * nt, 0.0); L.U0[k].assign(5 * nt, 0.0); L.W[k].assign(5 * nt, 0.0); L.Ux[k].assign(5 * nt, 0.0);
            L.g[k].assign(5 * n, 0.0); L.rate[k].assign(5 * n, 0.0); L.rr[k].assign(5 * n, 0.0);
        }
        L.dt.assign(n, 0.0);
        L.Btot.assign(3 * nt, 0.0); L.J.assign(3 * nt, 0.0); L.dpe.assign(3 * nt, 0.0); L.zero3.assign(3 * nt, 0.0);
        L.photo.assign(nt, 0.0); L.vst.assign(nt, 0.0);
    }
}

std::vector<long> Multigrid::global_sizes()
{
    std::vector<long> loc, glob;
    for (const MgLevel &L : lv_) loc.push_back(L.n_own);
    glob = loc;
    MPI_Allreduce(loc.data(), glob.data(), static_cast<int>(loc.size()), MPI_LONG, MPI_SUM, comm_);
    return glob;
}

void Multigrid::exchange(int l, std::initializer_list<std::pair<double *, int>> arrays)
{
    const MgLevel &L = lv_[l];
    const CellComplex &cx = lm_.mesh.cx;
    const Id nc = cx.n_cells(), no = cx.ent[3].n_owned;
    int M = 0;
    for (const auto &a : arrays) M += a.second;
    buf_.assign(static_cast<std::size_t>(nc) * M, 0.0);
    for (Id c = 0; c < no; ++c)
    {
        const Id cl = L.cell_cluster[c];
        if (cl < 0) continue;
        int off = 0;
        for (const auto &a : arrays)
        {
            for (int m = 0; m < a.second; ++m) buf_[static_cast<std::size_t>(c) * M + off + m] = a.first[static_cast<std::size_t>(cl) * a.second + m];
            off += a.second;
        }
    }
    ex_.forward(Dim::Cell, buf_.data(), M);
    for (Id j = 0; j < L.n_rem; ++j)
    {
        const Id g = L.rep_ghost[j];
        int off = 0;
        for (const auto &a : arrays)
        {
            for (int m = 0; m < a.second; ++m) a.first[static_cast<std::size_t>(L.n_own + j) * a.second + m] = buf_[static_cast<std::size_t>(g) * M + off + m];
            off += a.second;
        }
    }
}

void Multigrid::build_level(int l, const std::vector<Id> &own_id, Id n_own)
{
    MgLevel &L = lv_[l];
    const CellComplex &cx = lm_.mesh.cx;
    const Id nc = cx.n_cells(), no = cx.ent[3].n_owned;
    L.n_own = n_own;
    L.cell_cluster.assign(nc, -1);
    std::vector<double> idbuf(nc, -1.0);
    for (Id c = 0; c < no; ++c) idbuf[c] = static_cast<double>(own_id[c]);
    ex_.forward(Dim::Cell, idbuf.data(), 1);
    std::map<std::pair<int, long>, Id> rem;
    for (Id c = 0; c < nc; ++c)
    {
        if (c < no) { L.cell_cluster[c] = own_id[c]; continue; }
        if (idbuf[c] < 0) continue;
        const auto key = std::make_pair(cx.ent[3].owner[c], static_cast<long>(idbuf[c]));
        auto it = rem.find(key);
        if (it == rem.end())
        {
            it = rem.emplace(key, n_own + static_cast<Id>(rem.size())).first;
            L.rep_ghost.push_back(c);
        }
        L.cell_cluster[c] = it->second;
    }
    L.n_rem = static_cast<Id>(rem.size());
    L.vol.assign(L.n_tot(), 0.0);
    for (Id c = 0; c < no; ++c)
        if (L.cell_cluster[c] >= 0) L.vol[L.cell_cluster[c]] += geo_.cell_volume[c];
    exchange(l, {{L.vol.data(), 1}});

    struct TI { Id a, b; Vec3 A; };
    struct TB { Id a; int kind; Vec3 A; };
    std::vector<TI> inter;
    std::vector<TB> bnd;
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const FaceKind kind = bc_.face_kind[f];
        if (kind == FaceKind::Skip) continue;
        const auto k0 = cx.f2c.begin(f);
        const int nf = cx.f2c.count(f);
        if (kind == FaceKind::Interior && nf == 2)
        {
            const Id a = cx.f2c.idx[k0], b = cx.f2c.idx[k0 + 1];
            const Id ca = L.cell_cluster[a], cb = L.cell_cluster[b];
            if (ca < 0 || cb < 0 || ca == cb) continue;
            if (ca < n_own) inter.push_back({ca, cb, static_cast<double>(cx.f2c.sgn[k0]) * geo_.face_area[f]});
            if (cb < n_own) inter.push_back({cb, ca, static_cast<double>(cx.f2c.sgn[k0 + 1]) * geo_.face_area[f]});
        }
        else
            for (int q = 0; q < nf; ++q)
            {
                const Id c = cx.f2c.idx[k0 + q];
                if (cx.cell_material[c] != kFluid) continue;
                const Id cl = L.cell_cluster[c];
                if (cl >= 0 && cl < n_own) bnd.push_back({cl, static_cast<int>(kind), static_cast<double>(cx.f2c.sgn[k0 + q]) * geo_.face_area[f]});
            }
    }
    // Faces between the same two clusters (or of the same kind on the same cluster) are merged only when their normals agree (cosine > cos_tol):
    // a corner cluster must not get one "diagonal" boundary face of area |sum A| < sum |A|.  Closure (sum of all area vectors = 0) holds regardless.
    const double cos_tol = 0.9;
    auto merge_group = [&](std::vector<Vec3> &in, std::vector<Vec3> &out) {
        std::vector<Vec3> ref;
        for (const Vec3 &A : in)
        {
            const double a = norm(A);
            std::size_t hit = ref.size();
            for (std::size_t q = 0; q < ref.size(); ++q)
                if (a > 0 && dot(A, out[q]) > cos_tol * a * norm(out[q])) { hit = q; break; }
            if (hit == ref.size()) { ref.push_back(A); out.push_back(A); }
            else out[hit] += A;
        }
    };
    std::sort(inter.begin(), inter.end(), [](const TI &x, const TI &y) { return x.a != y.a ? x.a < y.a : x.b < y.b; });
    L.fptr.assign(n_own + 1, 0); L.fnb.clear(); L.fA.clear();
    for (std::size_t i = 0; i < inter.size();)
    {
        std::size_t j = i;
        std::vector<Vec3> in, out;
        for (; j < inter.size() && inter[j].a == inter[i].a && inter[j].b == inter[i].b; ++j) in.push_back(inter[j].A);
        merge_group(in, out);
        for (const Vec3 &A : out) { L.fnb.push_back(inter[i].b); L.fA.push_back(A); L.fptr[inter[i].a + 1]++; }
        i = j;
    }
    for (Id i = 0; i < n_own; ++i) L.fptr[i + 1] += L.fptr[i];
    std::sort(bnd.begin(), bnd.end(), [](const TB &x, const TB &y) { return x.a != y.a ? x.a < y.a : x.kind < y.kind; });
    L.bptr.assign(n_own + 1, 0); L.bkind.clear(); L.bA.clear();
    for (std::size_t i = 0; i < bnd.size();)
    {
        std::size_t j = i;
        std::vector<Vec3> in, out;
        for (; j < bnd.size() && bnd[j].a == bnd[i].a && bnd[j].kind == bnd[i].kind; ++j) in.push_back(bnd[j].A);
        merge_group(in, out);
        for (const Vec3 &A : out) { L.bkind.push_back(static_cast<FaceKind>(bnd[i].kind)); L.bA.push_back(A); L.bptr[bnd[i].a + 1]++; }
        i = j;
    }
    for (Id i = 0; i < n_own; ++i) L.bptr[i + 1] += L.bptr[i];
}

// heavy-face matching on the own clusters (cfg.mg_agglom_passes passes: pairs, then pairs of pairs, ...)
Id Multigrid::match(int l, std::vector<Id> &parent) const
{
    const MgLevel &L = lv_[l];
    const Id n = L.n_own;
    parent.resize(n);
    std::iota(parent.begin(), parent.end(), 0);
    if (n <= cfg_.mg_min_clusters) return n;
    std::vector<Id> grp(n);
    std::iota(grp.begin(), grp.end(), 0);
    Id ng = n;
    for (int pass = 0; pass < std::max(1, cfg_.mg_agglom_passes) && ng > cfg_.mg_min_clusters; ++pass)
    {
        std::vector<std::tuple<Id, Id, double>> ed;
        for (Id i = 0; i < n; ++i)
            for (auto k = L.fptr[i]; k < L.fptr[i + 1]; ++k)
            {
                const Id j = L.fnb[k];
                if (j >= n || grp[i] == grp[j]) continue;
                ed.emplace_back(grp[i], grp[j], norm(L.fA[k]));
            }
        std::sort(ed.begin(), ed.end(), [](const auto &x, const auto &y) { return std::get<0>(x) != std::get<0>(y) ? std::get<0>(x) < std::get<0>(y) : std::get<1>(x) < std::get<1>(y); });
        std::vector<std::int64_t> ptr(ng + 1, 0);
        std::vector<Id> nb;
        std::vector<double> w;
        for (std::size_t i = 0; i < ed.size();)
        {
            std::size_t j = i;
            double s = 0;
            for (; j < ed.size() && std::get<0>(ed[j]) == std::get<0>(ed[i]) && std::get<1>(ed[j]) == std::get<1>(ed[i]); ++j) s += std::get<2>(ed[j]);
            nb.push_back(std::get<1>(ed[i])); w.push_back(s); ptr[std::get<0>(ed[i]) + 1]++;
            i = j;
        }
        for (Id g = 0; g < ng; ++g) ptr[g + 1] += ptr[g];
        std::vector<Id> mate(ng, -1), rep(ng, -1);
        for (Id g = 0; g < ng; ++g)
        {
            if (mate[g] >= 0) continue;
            Id best = -1;
            double bw = 0;
            for (auto k = ptr[g]; k < ptr[g + 1]; ++k)
                if (mate[nb[k]] < 0 && w[k] > bw) { bw = w[k]; best = nb[k]; }
            if (best >= 0) { mate[g] = best; mate[best] = g; rep[g] = rep[best] = std::min(g, best); }
        }
        for (Id g = 0; g < ng; ++g) // unmatched groups join the neighbour they are most strongly connected to
        {
            if (rep[g] >= 0) continue;
            Id best = -1;
            double bw = 0;
            for (auto k = ptr[g]; k < ptr[g + 1]; ++k)
                if (rep[nb[k]] >= 0 && w[k] > bw) { bw = w[k]; best = nb[k]; }
            rep[g] = best >= 0 ? rep[best] : g;
        }
        std::vector<Id> label(ng, -1);
        Id nn = 0;
        for (Id g = 0; g < ng; ++g)
            if (rep[g] == g) label[g] = nn++;
        for (Id g = 0; g < ng; ++g) label[g] = label[rep[g]];
        if (nn >= ng) break;
        for (Id i = 0; i < n; ++i) grp[i] = label[grp[i]];
        ng = nn;
    }
    parent = grp;
    return ng;
}

void Multigrid::prim_all(int l, Id n, const double *const U[2], double *const W[2]) const
{
    (void)l;
    for (int k = 0; k < 2; ++k)
        for (Id i = 0; i < n; ++i)
        {
            const Prim w = to_prim(U[k] + 5 * static_cast<std::size_t>(i), cfg_.sc.gamma, cfg_.rho_floor, cfg_.p_floor);
            double *o = W[k] + 5 * static_cast<std::size_t>(i);
            o[0] = w.rho; o[1] = w.u; o[2] = w.v; o[3] = w.w; o[4] = w.p;
        }
}

namespace
{
inline Prim prim_at(const double *W, Id c) { return {W[5 * c], W[5 * c + 1], W[5 * c + 2], W[5 * c + 3], W[5 * c + 4]}; }
} // namespace

void Multigrid::flux_rhs(int l, const double *const W[2], double *const rhs[2]) const
{
    const MgLevel &L = lv_[l];
    for (int k = 0; k < 2; ++k)
        for (Id i = 0; i < L.n_own; ++i)
        {
            double acc[5] = {0, 0, 0, 0, 0};
            const Prim in = prim_at(W[k], i);
            for (auto q = L.fptr[i]; q < L.fptr[i + 1]; ++q)
            {
                const double A = norm(L.fA[q]);
                if (A <= 0) continue;
                const Vec3 n = (1.0 / A) * L.fA[q];
                const double nv[3] = {n.x, n.y, n.z};
                double F[5];
                rusanov(in, prim_at(W[k], L.fnb[q]), nv, cfg_.sc.gamma, F);
                for (int m = 0; m < 5; ++m) acc[m] += F[m] * A;
            }
            for (auto q = L.bptr[i]; q < L.bptr[i + 1]; ++q)
            {
                const double A = norm(L.bA[q]);
                if (A <= 0) continue;
                const Vec3 n = (1.0 / A) * L.bA[q];
                const double nv[3] = {n.x, n.y, n.z};
                double F[5];
                rusanov(in, ghost_state(L.bkind[q], k, in, n, fs_[k]), nv, cfg_.sc.gamma, F);
                for (int m = 0; m < 5; ++m) acc[m] += F[m] * A;
            }
            for (int m = 0; m < 5; ++m) rhs[k][5 * static_cast<std::size_t>(i) + m] = -acc[m] / L.vol[i];
        }
}

void Multigrid::dt_level(int l)
{
    MgLevel &L = lv_[l];
    for (Id i = 0; i < L.n_own; ++i)
    {
        const Vec3 b(L.Btot[3 * i], L.Btot[3 * i + 1], L.Btot[3 * i + 2]);
        const double B2 = dot(b, b);
        const double rho_tot = std::max(L.W[0][5 * static_cast<std::size_t>(i)] + L.W[1][5 * static_cast<std::size_t>(i)], cfg_.rho_floor_dt);
        double dt = 1e300;
        for (int s = 0; s < 2; ++s)
        {
            const double *w = L.W[s].data() + 5 * static_cast<std::size_t>(i);
            const double cs2 = cfg_.sc.gamma * w[4] / std::max(w[0], 1e-30);
            double denom = 0;
            auto add = [&](const Vec3 &Av) {
                const double A = norm(Av);
                if (A <= 0) return;
                const Vec3 n = (1.0 / A) * Av;
                const double vA2 = cfg_.sc.inver_MA2 * B2 / rho_tot, vAn2 = cfg_.sc.inver_MA2 * dot(b, n) * dot(b, n) / rho_tot;
                const double term = cs2 + vA2;
                const double cf = std::sqrt(0.5 * (term + std::sqrt(std::max(0.0, term * term - 4 * cs2 * vAn2))));
                denom += (std::abs(w[1] * n.x + w[2] * n.y + w[3] * n.z) + cf) * A;
            };
            for (auto q = L.fptr[i]; q < L.fptr[i + 1]; ++q) add(L.fA[q]);
            for (auto q = L.bptr[i]; q < L.bptr[i + 1]; ++q) add(L.bA[q]);
            if (denom > 0) dt = std::min(dt, cfg_.cfl * L.vol[i] / denom);
        }
        L.dt[i] = dt;
    }
}

void Multigrid::source_rhs(int l, double *const rhs[2])
{
    MgLevel &L = lv_[l];
    SourceFields sf;
    sf.photo = L.photo.data(); sf.vst = L.vst.data(); sf.dpe = L.dpe.data();
    for (int k = 0; k < 2; ++k)
        add_sources_n(L.n_own, nullptr, cfg_, k, L.W[0].data(), L.W[1].data(), L.zero3.data(), L.Btot.data(), L.J.data(), sf, 0.0, rhs[k], L.dt.data());
}

void Multigrid::floors(int l)
{
    MgLevel &L = lv_[l];
    for (int k = 0; k < 2; ++k)
        for (Id i = 0; i < L.n_own; ++i)
        {
            double *u = &L.U[k][5 * static_cast<std::size_t>(i)];
            const double ke = u[0] > 0 ? 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0] : 0.0;
            if (!(u[0] >= cfg_.rho_floor)) to_cons(Prim{cfg_.rho_floor, 0, 0, 0, cfg_.p_floor}, cfg_.sc.gamma, u);
            else if ((cfg_.sc.gamma - 1.0) * (u[4] - ke) < cfg_.p_floor) to_cons(to_prim(u, cfg_.sc.gamma, cfg_.rho_floor, cfg_.p_floor), cfg_.sc.gamma, u);
        }
}

void Multigrid::pseudo_step(int l, bool use_forcing)
{
    MgLevel &L = lv_[l];
    const Id n = L.n_own;
    double *Uk[2] = {L.U[0].data(), L.U[1].data()}, *Wk[2] = {L.W[0].data(), L.W[1].data()};
    prim_all(l, n, Uk, Wk);
    dt_level(l);
    const bool split = cfg_.sc.induce_coeff > 0 && cfg_.lorentz_force;
    if (split) gyrate_n(n, nullptr, cfg_, L.Btot.data(), 0.5, Uk[0], Uk[1], L.dt.data());
    exchange(l, {{Uk[0], 5}, {Uk[1], 5}});
    prim_all(l, L.n_tot(), Uk, Wk);
    std::vector<double> r0(5 * static_cast<std::size_t>(n)), r1(5 * static_cast<std::size_t>(n));
    double *rhs[2] = {r0.data(), r1.data()};
    flux_rhs(l, Wk, rhs);
    source_rhs(l, rhs);
    for (int k = 0; k < 2; ++k)
        for (Id i = 0; i < n; ++i)
            for (int m = 0; m < 5; ++m) Uk[k][5 * static_cast<std::size_t>(i) + m] += L.dt[i] * rhs[k][5 * static_cast<std::size_t>(i) + m];
    floors(l);
    if (split) gyrate_n(n, nullptr, cfg_, L.Btot.data(), 0.5, Uk[0], Uk[1], L.dt.data());
    if (use_forcing)
        for (int k = 0; k < 2; ++k)
            for (Id i = 0; i < n; ++i)
                for (int m = 0; m < 5; ++m) Uk[k][5 * static_cast<std::size_t>(i) + m] += L.dt[i] * L.g[k][5 * static_cast<std::size_t>(i) + m];
    if (cfg_.debug_res)
    {
        double dmin = 1e300, dmax = 0, umax = 0;
        long bad = 0;
        for (Id i = 0; i < n; ++i)
        {
            dmin = std::min(dmin, L.dt[i]); dmax = std::max(dmax, L.dt[i]);
            for (int k = 0; k < 2; ++k)
                for (int m = 0; m < 5; ++m)
                {
                    const double v = Uk[k][5 * static_cast<std::size_t>(i) + m];
                    if (!std::isfinite(v)) ++bad; else umax = std::max(umax, std::abs(v));
                }
        }
        if (bad || !std::isfinite(dmax) || dmax > 1e5) std::printf("[mg rank %d] level %d: nonfinite=%ld dt in [%.3e, %.3e] max|U|=%.3e forcing=%d\n", lm_.rank, l, bad, dmin, dmax, umax, int(use_forcing));
    }
}

void Multigrid::rate_of(int l, bool with_forcing)
{
    MgLevel &L = lv_[l];
    const Id n = L.n_own;
    for (int k = 0; k < 2; ++k) { L.Ux[k] = L.U[k]; std::swap(L.U[k], L.Ux[k]); } // Ux = original, U = scratch copy
    pseudo_step(l, false);
    for (int k = 0; k < 2; ++k)
    {
        for (Id i = 0; i < n; ++i)
            for (int m = 0; m < 5; ++m)
            {
                const std::size_t q = 5 * static_cast<std::size_t>(i) + m;
                L.rate[k][q] = (L.U[k][q] - L.Ux[k][q]) / L.dt[i] + (with_forcing ? L.g[k][q] : 0.0);
            }
        std::swap(L.U[k], L.Ux[k]);
    }
}

void Multigrid::restrict_fine(int l, const double *const Uf[2], double *const U[2]) const
{
    const MgLevel &L = lv_[l];
    const CellComplex &cx = lm_.mesh.cx;
    for (int k = 0; k < 2; ++k)
    {
        std::fill(U[k], U[k] + 5 * static_cast<std::size_t>(L.n_own), 0.0);
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            const Id cl = L.cell_cluster[c];
            if (cl < 0) continue;
            for (int m = 0; m < 5; ++m) U[k][5 * static_cast<std::size_t>(cl) + m] += geo_.cell_volume[c] * Uf[k][5 * static_cast<std::size_t>(c) + m];
        }
        for (Id i = 0; i < L.n_own; ++i)
            for (int m = 0; m < 5; ++m) U[k][5 * static_cast<std::size_t>(i) + m] /= L.vol[i];
    }
}

void Multigrid::restrict_up(int l)
{
    MgLevel &F = lv_[l], &C = lv_[l + 1];
    for (int k = 0; k < 2; ++k)
    {
        std::fill(C.U[k].begin(), C.U[k].begin() + 5 * static_cast<std::size_t>(C.n_own), 0.0);
        std::fill(C.rr[k].begin(), C.rr[k].end(), 0.0);
        for (Id i = 0; i < F.n_own; ++i)
        {
            const Id j = F.parent[i];
            for (int m = 0; m < 5; ++m)
            {
                C.U[k][5 * static_cast<std::size_t>(j) + m] += F.vol[i] * F.U[k][5 * static_cast<std::size_t>(i) + m];
                C.rr[k][5 * static_cast<std::size_t>(j) + m] += F.vol[i] * F.rate[k][5 * static_cast<std::size_t>(i) + m];
            }
        }
        for (Id j = 0; j < C.n_own; ++j)
            for (int m = 0; m < 5; ++m)
            {
                C.U[k][5 * static_cast<std::size_t>(j) + m] /= C.vol[j];
                C.rr[k][5 * static_cast<std::size_t>(j) + m] /= C.vol[j];
            }
    }
}

void Multigrid::prolong_down(int l, double relax)
{
    MgLevel &F = lv_[l - 1], &C = lv_[l];
    for (int k = 0; k < 2; ++k)
        for (Id i = 0; i < F.n_own; ++i)
        {
            const Id j = F.parent[i];
            double e[5];
            for (int m = 0; m < 5; ++m) e[m] = relax * (C.U[k][5 * static_cast<std::size_t>(j) + m] - C.U0[k][5 * static_cast<std::size_t>(j) + m]);
            double *u = &F.U[k][5 * static_cast<std::size_t>(i)];
            const double a = positivity_safeguard(u, e, cfg_.sc.gamma, 0.1, cfg_.mg_growth);
            for (int m = 0; m < 5; ++m) u[m] += a * e[m];
        }
}

void Multigrid::cycle(const double *const Uf[2], const double *const ratef[2], const double *Btot, const double *J, const double *dpe, const double *photo,
                      const double *vst, double *const dU[2])
{
    const double t0 = MPI_Wtime();
    const CellComplex &cx = lm_.mesh.cx;
    const int nl = n_levels();
    for (int k = 0; k < 2; ++k)
        for (Id c = 0; c < cx.n_cells(); ++c)
            for (int m = 0; m < 5; ++m) dU[k][5 * static_cast<std::size_t>(c) + m] = 0.0;
    if (nl < 2) return;
    // frozen fields on every level
    for (int l = 1; l < nl; ++l)
    {
        MgLevel &L = lv_[l];
        std::fill(L.Btot.begin(), L.Btot.end(), 0.0); std::fill(L.J.begin(), L.J.end(), 0.0); std::fill(L.dpe.begin(), L.dpe.end(), 0.0);
        std::fill(L.photo.begin(), L.photo.end(), 0.0); std::fill(L.vst.begin(), L.vst.end(), 0.0);
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            const Id cl = L.cell_cluster[c];
            if (cl < 0) continue;
            const double V = geo_.cell_volume[c];
            for (int a = 0; a < 3; ++a)
            {
                L.Btot[3 * cl + a] += V * Btot[3 * c + a]; L.J[3 * cl + a] += V * J[3 * c + a]; L.dpe[3 * cl + a] += V * dpe[3 * c + a];
            }
            L.photo[cl] += V * photo[c]; L.vst[cl] += V * vst[c];
        }
        for (Id i = 0; i < L.n_own; ++i)
        {
            for (int a = 0; a < 3; ++a) { L.Btot[3 * i + a] /= L.vol[i]; L.J[3 * i + a] /= L.vol[i]; L.dpe[3 * i + a] /= L.vol[i]; }
            L.photo[i] /= L.vol[i]; L.vst[i] /= L.vol[i];
        }
        exchange(l, {{L.Btot.data(), 3}, {L.J.data(), 3}, {L.dpe.data(), 3}, {L.photo.data(), 1}, {L.vst.data(), 1}});
    }
    // level 1 from the fine level
    {
        MgLevel &L = lv_[1];
        double *U[2] = {L.U[0].data(), L.U[1].data()}, *rr[2] = {L.rr[0].data(), L.rr[1].data()};
        restrict_fine(1, Uf, U);
        restrict_fine(1, ratef, rr);
    }
    // down sweep
    for (int l = 1; l < nl; ++l)
    {
        MgLevel &L = lv_[l];
        for (int k = 0; k < 2; ++k) L.U0[k] = L.U[k];
        rate_of(l, false);
        for (int k = 0; k < 2; ++k)
            for (std::size_t q = 0; q < L.g[k].size(); ++q) L.g[k][q] = L.rr[k][q] - L.rate[k][q];
        const int nsm = (l == nl - 1) ? cfg_.mg_coarse : cfg_.mg_pre;
        for (int s = 0; s < nsm; ++s) pseudo_step(l, true);
        if (l < nl - 1)
        {
            rate_of(l, true);
            restrict_up(l);
        }
    }
    // up sweep
    for (int l = nl - 1; l >= 2; --l)
    {
        prolong_down(l, 1.0);
        for (int s = 0; s < cfg_.mg_post; ++s) pseudo_step(l - 1, true);
    }
    // correction of level 1 -> fine cells
    {
        const MgLevel &L = lv_[1];
        for (int k = 0; k < 2; ++k)
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            {
                const Id cl = L.cell_cluster[c];
                if (cl < 0) continue;
                for (int m = 0; m < 5; ++m)
                    dU[k][5 * static_cast<std::size_t>(c) + m] = cfg_.mg_relax * (L.U[k][5 * static_cast<std::size_t>(cl) + m] - L.U0[k][5 * static_cast<std::size_t>(cl) + m]);
            }
    }
    t_cycle_ += MPI_Wtime() - t0;
}
} // namespace mercury
