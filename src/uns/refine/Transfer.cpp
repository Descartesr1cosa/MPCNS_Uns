// Coarse -> fine transfer of the face flux B (divergence free per fine cell) and of cell fields.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include <unordered_map>

#include "uns/refine/Refine.h"

namespace uns
{
void detail_lsq_gradient(const CellComplex &cx, const Geometry &g, const std::vector<double> &u, int ncomp, int comp,
                         std::vector<Vec3> &grad, std::vector<double> *umin, std::vector<double> *umax);

namespace
{
// area-weighted least-squares cell vector b with  b . n_f = flux_f / A_f  (outward n_f)
std::vector<Vec3> reconstruct_cell_vector(const Mesh &m, const Geometry &g, const std::vector<double> &flux)
{
    const CellComplex &cx = m.cx;
    std::vector<Vec3> b(cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        double M[3][3] = {{0}}, r[3] = {0};
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            const double A = g.face_area_mag[f];
            if (A <= 0) continue;
            const Vec3 n = (cx.c2f.sgn[k] / A) * g.face_area[f];
            const double nv[3] = {n.x, n.y, n.z};
            const double phi = cx.c2f.sgn[k] * flux[f];
            for (int i = 0; i < 3; ++i)
            {
                r[i] += nv[i] * phi; // sum_f A_f n (flux/A)  with weight A_f  ->  n * flux
                for (int j = 0; j < 3; ++j) M[i][j] += A * nv[i] * nv[j];
            }
        }
        const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        const double sc = std::max({M[0][0], M[1][1], M[2][2], 1e-300});
        if (std::abs(det) < 1e-12 * sc * sc * sc) continue;
        double inv[3][3];
        inv[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det;
        inv[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / det;
        inv[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
        inv[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) / det;
        inv[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
        inv[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / det;
        inv[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det;
        inv[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / det;
        inv[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
        b[c] = {inv[0][0] * r[0] + inv[0][1] * r[1] + inv[0][2] * r[2], inv[1][0] * r[0] + inv[1][1] * r[1] + inv[1][2] * r[2],
                inv[2][0] * r[0] + inv[2][1] * r[1] + inv[2][2] * r[2]};
    }
    return b;
}

int sign_of(const Incidence &t, Id row, Id col)
{
    for (auto k = t.begin(row); k < t.end(row); ++k)
        if (t.idx[k] == col) return t.sgn[k];
    return 0;
}

// Make every cell of one coarse region divergence free by changing interior-face fluxes only:
// minimum-norm correction (CG on the cell graph Laplacian) followed by a spanning-tree sweep.
void project_region(const CellComplex &fx, const Geometry &fg, const std::vector<Id> &interior_mask_face_coarse, Id c0, Id c1,
                    std::vector<double> &phi)
{
    const Id nk = c1 - c0;
    if (nk <= 1) return;
    std::unordered_map<Id, int> lid; // interior fine face -> local index
    std::vector<Id> faces;
    for (Id c = c0; c < c1; ++c)
        for (auto k = fx.c2f.begin(c); k < fx.c2f.end(c); ++k)
        {
            const Id f = fx.c2f.idx[k];
            if (interior_mask_face_coarse[f] == -1 && lid.emplace(f, static_cast<int>(faces.size())).second) faces.push_back(f);
        }
    const int nf = static_cast<int>(faces.size());
    struct Ends { int c[2]; int s[2]; };
    std::vector<Ends> ends(nf);
    for (int i = 0; i < nf; ++i)
    {
        const Id f = faces[i];
        for (int q = 0; q < 2; ++q)
        {
            const auto k = fx.f2c.begin(f) + q;
            ends[i].c[q] = fx.f2c.idx[k] - c0;
            ends[i].s[q] = fx.f2c.sgn[k];
        }
    }
    auto net_out = [&](Id c) {
        double s = 0;
        for (auto k = fx.c2f.begin(c); k < fx.c2f.end(c); ++k) s += fx.c2f.sgn[k] * phi[fx.c2f.idx[k]];
        return s;
    };
    std::vector<double> rhs(nk), x(nk, 0), res(nk), p(nk), Ap(nk), diag(nk, 0), W(nf);
    double scale = 0;
    for (Id c = c0; c < c1; ++c)
    {
        rhs[c - c0] = -net_out(c);
        for (auto k = fx.c2f.begin(c); k < fx.c2f.end(c); ++k) scale = std::max(scale, std::abs(phi[fx.c2f.idx[k]]));
    }
    for (int i = 0; i < nf; ++i)
    {
        W[i] = fg.face_area_mag[faces[i]];
        diag[ends[i].c[0]] += W[i];
        diag[ends[i].c[1]] += W[i];
    }
    const double mean = std::accumulate(rhs.begin(), rhs.end(), 0.0) / nk;
    for (double &r : rhs) r -= mean; // compatibility (coarse d2 B = 0)
    auto matvec = [&](const std::vector<double> &v, std::vector<double> &y) {
        std::fill(y.begin(), y.end(), 0.0);
        for (int i = 0; i < nf; ++i)
        {
            const double t = W[i] * (ends[i].s[0] * v[ends[i].c[0]] + ends[i].s[1] * v[ends[i].c[1]]);
            y[ends[i].c[0]] += ends[i].s[0] * t;
            y[ends[i].c[1]] += ends[i].s[1] * t;
        }
    };
    res = rhs;
    std::vector<double> z(nk);
    for (int i = 0; i < nk; ++i) z[i] = diag[i] > 0 ? res[i] / diag[i] : 0;
    p = z;
    double rz = 0;
    for (int i = 0; i < nk; ++i) rz += res[i] * z[i];
    const int maxit = 20 * nk + 200;
    for (int it = 0; it < maxit && rz > 0; ++it)
    {
        double rmax = 0;
        for (double r : res) rmax = std::max(rmax, std::abs(r));
        if (rmax <= 1e-15 * scale) break;
        matvec(p, Ap);
        double pAp = 0;
        for (int i = 0; i < nk; ++i) pAp += p[i] * Ap[i];
        if (pAp <= 0) break;
        const double al = rz / pAp;
        for (int i = 0; i < nk; ++i) { x[i] += al * p[i]; res[i] -= al * Ap[i]; }
        for (int i = 0; i < nk; ++i) z[i] = diag[i] > 0 ? res[i] / diag[i] : 0;
        double rz2 = 0;
        for (int i = 0; i < nk; ++i) rz2 += res[i] * z[i];
        const double be = rz2 / rz;
        rz = rz2;
        for (int i = 0; i < nk; ++i) p[i] = z[i] + be * p[i];
    }
    for (int i = 0; i < nf; ++i) phi[faces[i]] += W[i] * (ends[i].s[0] * x[ends[i].c[0]] + ends[i].s[1] * x[ends[i].c[1]]);

    // spanning-tree sweep: removes the remaining (round-off sized) residuals
    std::vector<std::vector<int>> adj(nk);
    for (int i = 0; i < nf; ++i) { adj[ends[i].c[0]].push_back(i); adj[ends[i].c[1]].push_back(i); }
    std::vector<int> parent_face(nk, -1), order;
    std::vector<char> seen(nk, 0);
    std::queue<int> q;
    q.push(0); seen[0] = 1;
    while (!q.empty())
    {
        const int c = q.front(); q.pop();
        order.push_back(c);
        for (int fi : adj[c])
        {
            const int o = ends[fi].c[0] == c ? ends[fi].c[1] : ends[fi].c[0];
            if (!seen[o]) { seen[o] = 1; parent_face[o] = fi; q.push(o); }
        }
    }
    std::vector<double> out(nk);
    for (Id c = c0; c < c1; ++c) out[c - c0] = net_out(c);
    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const int c = *it;
        const int fi = parent_face[c];
        if (fi < 0) continue;
        const int side = ends[fi].c[0] == c ? 0 : 1;
        const double d = -out[c] * ends[fi].s[side]; // flux change on the face that zeroes cell c
        phi[faces[fi]] += d;
        out[c] = 0;
        out[ends[fi].c[1 - side]] += ends[fi].s[1 - side] * d;
    }
}
} // namespace

std::vector<double> transfer_B(const Mesh &coarse, const Geometry &cg, const std::vector<double> &coarse_flux, const RefineResult &rr,
                               const Geometry &fg, TransferReport *rep)
{
    const CellComplex &cx = coarse.cx, &fx = rr.fine.cx;
    const std::vector<Vec3> b = reconstruct_cell_vector(coarse, cg, coarse_flux);
    double flux_scale = 1e-300;
    for (double v : coarse_flux) flux_scale = std::max(flux_scale, std::abs(v));
    std::vector<double> phi(fx.n_faces(), 0.0);

    // faces lying on coarse faces: group the tiles per coarse face
    std::vector<std::vector<Id>> tiles(cx.n_faces());
    for (Id f = 0; f < fx.n_faces(); ++f)
        if (rr.face_coarse_face[f] >= 0) tiles[rr.face_coarse_face[f]].push_back(f);
    auto sigma = [&](Id f, Id F) {
        const Id L = rr.face_owner_cell[f];
        return sign_of(fx.c2f, L, f) * sign_of(cx.c2f, rr.parent_cell[L], F);
    };
    for (Id F = 0; F < cx.n_faces(); ++F)
    {
        if (tiles[F].empty()) continue;
        Vec3 bf;
        int nb = 0;
        for (auto k = cx.f2c.begin(F); k < cx.f2c.end(F); ++k) { bf += b[cx.f2c.idx[k]]; ++nb; }
        bf *= 1.0 / nb;
        double raw_sum = 0, area_sum = 0;
        std::vector<double> raw(tiles[F].size());
        for (std::size_t i = 0; i < tiles[F].size(); ++i)
        {
            const Id f = tiles[F][i];
            raw[i] = dot(bf, sigma(f, F) * fg.face_area[f]); // flux in the coarse canonical direction
            raw_sum += raw[i];
            area_sum += fg.face_area_mag[f];
        }
        const double defect = coarse_flux[F] - raw_sum; // spread by area: tile fluxes sum to the coarse flux exactly
        double check = 0;
        for (std::size_t i = 0; i < tiles[F].size(); ++i)
        {
            const Id f = tiles[F][i];
            const double v = raw[i] + defect * fg.face_area_mag[f] / area_sum;
            phi[f] = sigma(f, F) * v;
            check += v;
        }
        if (rep) rep->max_tile_defect = std::max(rep->max_tile_defect, std::abs(check - coarse_flux[F]) / flux_scale);
    }
    // interior faces of each coarse cell: smooth initial guess, then divergence-free projection
    for (Id f = 0; f < fx.n_faces(); ++f)
        if (rr.face_coarse_face[f] < 0) phi[f] = dot(b[rr.parent_cell[rr.face_owner_cell[f]]], fg.face_area[f]);
    for (Id c0 = 0; c0 < fx.n_cells();)
    {
        Id c1 = c0;
        while (c1 < fx.n_cells() && rr.parent_cell[c1] == rr.parent_cell[c0]) ++c1;
        project_region(fx, fg, rr.face_coarse_face, c0, c1, phi);
        c0 = c1;
    }
    if (rep)
        for (Id c = 0; c < fx.n_cells(); ++c)
        {
            double s = 0, a = 0;
            for (auto k = fx.c2f.begin(c); k < fx.c2f.end(c); ++k) { s += fx.c2f.sgn[k] * phi[fx.c2f.idx[k]]; a += std::abs(phi[fx.c2f.idx[k]]); }
            rep->max_rel_div = std::max(rep->max_rel_div, a > 0 ? std::abs(s) / a : 0.0);
        }
    return phi;
}

std::vector<double> prolong_cells(const Mesh &coarse, const Geometry &cg, const std::vector<double> &U, int ncomp, const RefineResult &rr,
                                  const Geometry &fg, bool limit)
{
    const Id nf = rr.fine.cx.n_cells();
    std::vector<double> out(static_cast<std::size_t>(nf) * ncomp);
    std::vector<Vec3> grad;
    std::vector<double> umin, umax;
    // children ranges per coarse cell
    std::vector<Id> first(coarse.cx.n_cells() + 1, 0);
    for (Id c = 0; c < nf; ++c) ++first[rr.parent_cell[c] + 1];
    for (Id c = 0; c < coarse.cx.n_cells(); ++c) first[c + 1] += first[c];
    for (int comp = 0; comp < ncomp; ++comp)
    {
        detail_lsq_gradient(coarse.cx, cg, U, ncomp, comp, grad, &umin, &umax);
        for (Id C = 0; C < coarse.cx.n_cells(); ++C)
        {
            const double uc = U[static_cast<std::size_t>(C) * ncomp + comp];
            const Id a = first[C], b = first[C + 1];
            if (b - a == 1) { out[static_cast<std::size_t>(a) * ncomp + comp] = uc; continue; }
            double phi = 1.0;
            if (limit)
                for (Id f = a; f < b; ++f)
                {
                    const double d = dot(grad[C], fg.cell_center[f] - cg.cell_center[C]);
                    if (d > 1e-300) phi = std::min(phi, (umax[C] - uc) / d);
                    else if (d < -1e-300) phi = std::min(phi, (umin[C] - uc) / d);
                }
            phi = std::max(0.0, std::min(1.0, phi));
            double vsum = 0, isum = 0;
            for (Id f = a; f < b; ++f)
            {
                const double inc = phi * dot(grad[C], fg.cell_center[f] - cg.cell_center[C]);
                out[static_cast<std::size_t>(f) * ncomp + comp] = inc;
                vsum += fg.cell_volume[f];
                isum += fg.cell_volume[f] * inc;
            }
            const double mean = isum / vsum; // conservation: the volume-weighted increment vanishes
            for (Id f = a; f < b; ++f) out[static_cast<std::size_t>(f) * ncomp + comp] += uc - mean;
        }
    }
    return out;
}
} // namespace uns
