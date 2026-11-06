#include "fluid/Reconstruction.h"

#include <algorithm>

namespace mhd
{
using namespace uns;

Prim ghost_state(FaceKind kind, int species, const Prim &in, const Vec3 &n, const Prim &fs)
{
    switch (kind)
    {
    case FaceKind::Wall:
    {
        if (species == 1) return Prim{1e-12, 0, 0, 0, 1e-12}; // absorbing
        const double un = in.u * n.x + in.v * n.y + in.w * n.z;
        Prim g = in;
        g.u -= 2 * un * n.x; g.v -= 2 * un * n.y; g.w -= 2 * un * n.z;
        return g;
    }
    case FaceKind::Inflow: return fs;
    default: return in;
    }
}

Stencil::Stencil(const LocalMesh &lm, const Geometry &geo, const BoundaryInfo &bc)
{
    const CellComplex &cx = lm.mesh.cx;
    const Id no = cx.ent[3].n_owned;
    ptr_.assign(no + 1, 0);
    inv_.assign(9 * static_cast<std::size_t>(no), 0.0);
    valid_.assign(no, 0);
    for (Id c = 0; c < no; ++c)
    {
        if (cx.cell_material[c] == kSolid) { ptr_[c + 1] = ptr_[c]; continue; }
        double M[3][3] = {{0}};
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            Vec3 d;
            Id nb;
            const FaceKind kind = bc.face_kind[f];
            if (kind == FaceKind::Skip) continue;
            if (cx.f2c.count(f) == 2 && kind == FaceKind::Interior)
            {
                nb = cx.f2c.idx[cx.f2c.begin(f)] == c ? cx.f2c.idx[cx.f2c.begin(f) + 1] : cx.f2c.idx[cx.f2c.begin(f)];
                d = geo.cell_center[nb] - geo.cell_center[c];
            }
            else // boundary / interface face: mirrored virtual neighbour
            {
                const Vec3 n = (cx.c2f.sgn[k] / geo.face_area_mag[f]) * geo.face_area[f];
                d = (2.0 * dot(geo.face_center[f] - geo.cell_center[c], n)) * n;
                nb = -(f + 1);
            }
            const double w = 1.0 / dot(d, d);
            nb_.push_back(nb); w_.push_back(w); d_.push_back(d);
            const double dv[3] = {d.x, d.y, d.z};
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) M[a][b] += w * dv[a] * dv[b];
        }
        ptr_[c + 1] = static_cast<std::int64_t>(nb_.size());
        const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        const double sc = std::max({M[0][0], M[1][1], M[2][2], 1e-300});
        if (std::abs(det) < 1e-9 * sc * sc * sc) continue;
        double* I = &inv_[9 * static_cast<std::size_t>(c)];
        I[0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det;
        I[1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / det;
        I[2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
        I[3] = I[1];
        I[4] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
        I[5] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / det;
        I[6] = I[2];
        I[7] = I[5];
        I[8] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
        valid_[c] = 1;
    }
}

void Stencil::limited_gradients(const LocalMesh &lm, const Geometry &geo, const BoundaryInfo &bc, const double *W,
                                int species, const Prim &fs, bool limit, double *grad, const GhostFn *ghost) const
{
    const CellComplex &cx = lm.mesh.cx;
    const Id no = cx.ent[3].n_owned;
    std::fill(grad, grad + 15 * static_cast<std::size_t>(cx.n_cells()), 0.0);
    auto prim_at = [&](Id c) { return Prim{W[5 * c], W[5 * c + 1], W[5 * c + 2], W[5 * c + 3], W[5 * c + 4]}; };
    for (Id c = 0; c < no; ++c)
    {
        if (!valid_[c]) continue;
        const Prim wc = prim_at(c);
        const double wcv[5] = {wc.rho, wc.u, wc.v, wc.w, wc.p};
        double r[5][3] = {{0}}, lo[5], hi[5];
        for (int m = 0; m < 5; ++m) lo[m] = hi[m] = wcv[m];
        for (auto k = ptr_[c]; k < ptr_[c + 1]; ++k)
        {
            Prim wn;
            if (nb_[k] >= 0) wn = prim_at(nb_[k]);
            else
            {
                const Id f = -(nb_[k] + 1);
                // outward normal of this cell through face f
                std::int64_t kk = cx.c2f.begin(c);
                while (cx.c2f.idx[kk] != f) ++kk;
                const Vec3 n = (cx.c2f.sgn[kk] / geo.face_area_mag[f]) * geo.face_area[f];
                wn = (ghost && *ghost) ? (*ghost)(f, species, wc, n) : ghost_state(bc.face_kind[f], species, wc, n, fs);
            }
            const double wnv[5] = {wn.rho, wn.u, wn.v, wn.w, wn.p};
            for (int m = 0; m < 5; ++m)
            {
                const double dq = wnv[m] - wcv[m];
                r[m][0] += w_[k] * d_[k].x * dq; r[m][1] += w_[k] * d_[k].y * dq; r[m][2] += w_[k] * d_[k].z * dq;
                lo[m] = std::min(lo[m], wnv[m]);
                hi[m] = std::max(hi[m], wnv[m]);
            }
        }
        const double *I = &inv_[9 * static_cast<std::size_t>(c)];
        for (int m = 0; m < 5; ++m)
        {
            double g[3] = {I[0] * r[m][0] + I[1] * r[m][1] + I[2] * r[m][2], I[3] * r[m][0] + I[4] * r[m][1] + I[5] * r[m][2],
                           I[6] * r[m][0] + I[7] * r[m][1] + I[8] * r[m][2]};
            double psi = 1.0;
            if (limit)
                for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
                {
                    const Vec3 dx = geo.face_center[cx.c2f.idx[kf]] - geo.cell_center[c];
                    const double dl = g[0] * dx.x + g[1] * dx.y + g[2] * dx.z;
                    if (dl > 1e-300) psi = std::min(psi, (hi[m] - wcv[m]) / dl);
                    else if (dl < -1e-300) psi = std::min(psi, (lo[m] - wcv[m]) / dl);
                }
            psi = std::max(0.0, std::min(1.0, psi));
            double *out = grad + 15 * static_cast<std::size_t>(c) + 3 * m;
            out[0] = psi * g[0]; out[1] = psi * g[1]; out[2] = psi * g[2];
        }
    }
}
} // namespace mhd
