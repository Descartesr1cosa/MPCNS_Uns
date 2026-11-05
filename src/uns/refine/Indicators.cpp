#include <algorithm>
#include <cmath>

#include "uns/refine/Refine.h"

namespace uns
{
namespace
{
// solve the 3x3 system M g = r (M symmetric); returns false if singular
bool solve3(const double M[3][3], const double r[3], double g[3])
{
    const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                       M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    double scale = std::max({std::abs(M[0][0]), std::abs(M[1][1]), std::abs(M[2][2]), 1e-300});
    if (std::abs(det) < 1e-12 * scale * scale * scale) return false;
    g[0] = (r[0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (r[1] * M[2][2] - M[1][2] * r[2]) +
            M[0][2] * (r[1] * M[2][1] - M[1][1] * r[2])) / det;
    g[1] = (M[0][0] * (r[1] * M[2][2] - M[1][2] * r[2]) - r[0] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
            M[0][2] * (M[1][0] * r[2] - r[1] * M[2][0])) / det;
    g[2] = (M[0][0] * (M[1][1] * r[2] - r[1] * M[2][1]) - M[0][1] * (M[1][0] * r[2] - r[1] * M[2][0]) +
            r[0] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
    return true;
}
} // namespace

// Inverse-distance weighted least-squares gradient over the face neighbours.
void detail_lsq_gradient(const CellComplex &cx, const Geometry &g, const std::vector<double> &u, int ncomp, int comp,
                         std::vector<Vec3> &grad, std::vector<double> *umin, std::vector<double> *umax)
{
    const Id nc = cx.n_cells();
    grad.assign(nc, Vec3());
    if (umin) umin->assign(nc, 0);
    if (umax) umax->assign(nc, 0);
    for (Id c = 0; c < nc; ++c)
    {
        double M[3][3] = {{0}}, r[3] = {0};
        const double uc = u[static_cast<std::size_t>(c) * ncomp + comp];
        double lo = uc, hi = uc;
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            const Id f = cx.c2f.idx[kf];
            for (auto q = cx.f2c.begin(f); q < cx.f2c.end(f); ++q)
            {
                const Id n = cx.f2c.idx[q];
                if (n == c) continue;
                const Vec3 d = g.cell_center[n] - g.cell_center[c];
                const double w = 1.0 / dot(d, d);
                const double dv[3] = {d.x, d.y, d.z};
                const double un = u[static_cast<std::size_t>(n) * ncomp + comp];
                lo = std::min(lo, un);
                hi = std::max(hi, un);
                for (int i = 0; i < 3; ++i)
                {
                    r[i] += w * dv[i] * (un - uc);
                    for (int j = 0; j < 3; ++j) M[i][j] += w * dv[i] * dv[j];
                }
            }
        }
        double gr[3];
        if (solve3(M, r, gr)) grad[c] = {gr[0], gr[1], gr[2]};
        if (umin) (*umin)[c] = lo;
        if (umax) (*umax)[c] = hi;
    }
}

std::vector<int> compute_target_levels(const Mesh &coarse, const Geometry &g, const RefineSpec &spec, const CoarseSolution *sol)
{
    const CellComplex &cx = coarse.cx;
    const Id nc = cx.n_cells();
    std::vector<int> lev(nc, 0);
    if (sol)
    {
        // cell size h = 6V / sum|A|
        std::vector<double> h(nc);
        for (Id c = 0; c < nc; ++c)
        {
            double a = 0;
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k) a += g.face_area_mag[cx.c2f.idx[k]];
            h[c] = 6.0 * g.cell_volume[c] / std::max(a, 1e-300);
        }
        std::vector<Vec3> grad;
        for (const IndicatorRule &rule : spec.indicators)
        {
            std::vector<double> ind(nc, 0.0);
            if (rule.kind == IndicatorRule::Kind::Current)
            {
                if (sol->J_cell.size() != static_cast<std::size_t>(nc) || sol->B_cell.size() != static_cast<std::size_t>(nc)) continue;
                for (Id c = 0; c < nc; ++c) ind[c] = h[c] * norm(sol->J_cell[c]) / std::max(norm(sol->B_cell[c]), 1e-300);
            }
            else
            {
                const std::vector<double> &f = rule.kind == IndicatorRule::Kind::Density ? sol->rho : sol->p;
                if (f.size() != static_cast<std::size_t>(nc)) continue;
                detail_lsq_gradient(cx, g, f, 1, 0, grad, nullptr, nullptr);
                for (Id c = 0; c < nc; ++c) ind[c] = h[c] * norm(grad[c]) / std::max(std::abs(f[c]), 1e-300);
            }
            for (Id c = 0; c < nc; ++c)
                if (ind[c] > rule.threshold) lev[c] = std::max(lev[c], rule.level);
        }
    }
    for (Id c = 0; c < nc; ++c)
    {
        lev[c] = std::min(lev[c], spec.max_level);
        if (cx.cell_material[c] == kSolid) lev[c] = std::min(lev[c], spec.solid_max_level);
    }
    return lev;
}
} // namespace uns
