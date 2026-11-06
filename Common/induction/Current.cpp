#include "induction/Current.h"

#include <algorithm>

namespace mhd
{
using namespace uns;

void cell_current(const LocalMesh &lm, const Geometry &geo, const double *B, double *J)
{
    const CellComplex &cx = lm.mesh.cx;
    std::fill(J, J + 3 * static_cast<std::size_t>(cx.n_cells()), 0.0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        Vec3 sum;
        const Vec3 bc(B[3 * c], B[3 * c + 1], B[3 * c + 2]);
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            Vec3 bf = bc;
            if (cx.f2c.count(f) == 2)
            {
                const Id o = cx.f2c.idx[cx.f2c.begin(f)] == c ? cx.f2c.idx[cx.f2c.begin(f) + 1] : cx.f2c.idx[cx.f2c.begin(f)];
                bf = 0.5 * (bc + Vec3(B[3 * o], B[3 * o + 1], B[3 * o + 2]));
            }
            sum += cross(static_cast<double>(cx.c2f.sgn[k]) * geo.face_area[f], bf);
        }
        const double invV = 1.0 / geo.cell_volume[c];
        J[3 * c] = sum.x * invV; J[3 * c + 1] = sum.y * invV; J[3 * c + 2] = sum.z * invV;
    }
}
} // namespace mhd
