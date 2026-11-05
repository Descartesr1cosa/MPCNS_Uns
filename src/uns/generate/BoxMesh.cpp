#include "uns/generate/HexMesh.h"

namespace uns
{
HexMesh make_box(const BoxSpec &s)
{
    HexMesh h;
    const int nx = s.nx, ny = s.ny, nz = s.nz;
    auto nid = [&](int i, int j, int k) { return static_cast<Id>((k * (ny + 1) + j) * (nx + 1) + i); };
    h.x.resize(static_cast<std::size_t>(nx + 1) * (ny + 1) * (nz + 1));
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i)
                h.x[nid(i, j, k)] = {s.lo.x + (s.hi.x - s.lo.x) * i / nx, s.lo.y + (s.hi.y - s.lo.y) * j / ny,
                                     s.lo.z + (s.hi.z - s.lo.z) * k / nz};
    h.tag_names = {"", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                std::array<Id, 8> v;
                for (int q = 0; q < 8; ++q) v[q] = nid(i + (q & 1), j + ((q >> 1) & 1), k + ((q >> 2) & 1));
                h.hex.push_back(v);
                h.material.push_back(kFluid);
                if (i == 0) h.tagged.push_back({{v[0], v[2], v[4], v[6]}, 1});
                if (i == nx - 1) h.tagged.push_back({{v[1], v[3], v[5], v[7]}, 2});
                if (j == 0) h.tagged.push_back({{v[0], v[1], v[4], v[5]}, 3});
                if (j == ny - 1) h.tagged.push_back({{v[2], v[3], v[6], v[7]}, 4});
                if (k == 0) h.tagged.push_back({{v[0], v[1], v[2], v[3]}, 5});
                if (k == nz - 1) h.tagged.push_back({{v[4], v[5], v[6], v[7]}, 6});
            }
    return h;
}
} // namespace uns
