#include "uns/geometry/Reconstruct.h"

#include <algorithm>
#include <cmath>

namespace uns
{
std::vector<Vec3> reconstruct_cell_vector(const CellComplex &cx, const Geometry &g, const double *flux, Id ncells)
{
    std::vector<Vec3> b(ncells);
    for (Id c = 0; c < ncells; ++c)
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

std::vector<Vec3> reconstruct_cell_vector(const Mesh &m, const Geometry &g, const std::vector<double> &flux)
{
    return reconstruct_cell_vector(m.cx, g, flux.data(), m.cx.n_cells());
}
} // namespace uns
