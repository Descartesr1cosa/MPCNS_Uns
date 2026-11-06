#include "fluid/FluxRHS.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mhd
{
using namespace uns;

namespace
{
inline Prim prim_at(const double *W, Id c) { return {W[5 * c], W[5 * c + 1], W[5 * c + 2], W[5 * c + 3], W[5 * c + 4]}; }

// linear reconstruction of cell c at point x, with positivity fallback to the cell value
inline Prim face_state(const double *W, const double *grad, Id c, const Vec3 &dx)
{
    Prim w = prim_at(W, c);
    const double *g = grad + 15 * static_cast<std::size_t>(c);
    Prim r;
    r.rho = w.rho + g[0] * dx.x + g[1] * dx.y + g[2] * dx.z;
    r.u = w.u + g[3] * dx.x + g[4] * dx.y + g[5] * dx.z;
    r.v = w.v + g[6] * dx.x + g[7] * dx.y + g[8] * dx.z;
    r.w = w.w + g[9] * dx.x + g[10] * dx.y + g[11] * dx.z;
    r.p = w.p + g[12] * dx.x + g[13] * dx.y + g[14] * dx.z;
    return (r.rho > 0 && r.p > 0) ? r : w;
}
} // namespace

void fluid_flux_rhs(const LocalMesh &lm, const Geometry &geo, const BoundaryInfo &bc, const Config &cfg, int species, const Prim &fs, const double *W,
                    const double *grad, double *rhs, const GhostFn *ghost)
{
    const CellComplex &cx = lm.mesh.cx;
    const Id no = cx.ent[3].n_owned;
    std::fill(rhs, rhs + 5 * static_cast<std::size_t>(cx.n_cells()), 0.0);
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const FaceKind kind = bc.face_kind[f];
        if (kind == FaceKind::Skip) continue;
        const auto k0 = cx.f2c.begin(f);
        const int nc = cx.f2c.count(f);
        // owned fluid cells touching the face
        bool any_owned = false;
        for (int q = 0; q < nc; ++q) any_owned |= cx.f2c.idx[k0 + q] < no && cx.cell_material[cx.f2c.idx[k0 + q]] == kFluid;
        if (!any_owned) continue;
        const double A = geo.face_area_mag[f];
        const Vec3 n = (1.0 / A) * geo.face_area[f]; // canonical normal, L -> R
        const double nv[3] = {n.x, n.y, n.z};
        Prim L, R;
        Id cL = -1, cR = -1;
        if (nc == 2 && kind == FaceKind::Interior)
        {
            cL = cx.f2c.idx[k0]; cR = cx.f2c.idx[k0 + 1];
            if (cx.f2c.sgn[k0] < 0) std::swap(cL, cR);
            L = face_state(W, grad, cL, geo.face_center[f] - geo.cell_center[cL]);
            R = face_state(W, grad, cR, geo.face_center[f] - geo.cell_center[cR]);
        }
        else
        {
            // one fluid cell + a boundary / interface (wall) ghost
            Id c = cx.f2c.idx[k0];
            int sgn = cx.f2c.sgn[k0];
            if (nc == 2 && cx.cell_material[c] == kSolid) { c = cx.f2c.idx[k0 + 1]; sgn = cx.f2c.sgn[k0 + 1]; }
            const Prim in = face_state(W, grad, c, geo.face_center[f] - geo.cell_center[c]);
            const Prim gh = (ghost && *ghost) ? (*ghost)(f, species, in, (double(sgn)) * n) : ghost_state(kind, species, in, (double(sgn)) * n, fs);
            if (sgn > 0) { cL = c; L = in; R = gh; } else { cR = c; R = in; L = gh; }
        }
        double F[5];
        rusanov(L, R, nv, cfg.sc.gamma, F);
        for (int m = 0; m < 5; ++m) F[m] *= A;
        if (cL >= 0 && cL < no) for (int m = 0; m < 5; ++m) rhs[5 * static_cast<std::size_t>(cL) + m] -= F[m] / geo.cell_volume[cL];
        if (cR >= 0 && cR < no) for (int m = 0; m < 5; ++m) rhs[5 * static_cast<std::size_t>(cR) + m] += F[m] / geo.cell_volume[cR];
    }
}
} // namespace mhd
