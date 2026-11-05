// Cubed-sphere shell.  Nodes are keyed by integer cube-surface lattice points
// (a,b,c) in [-N,N]^3 with max|.| = N and a radial layer index, so the panel
// seams and the 8 cube corners (valence 3) merge without any coordinate matching.
#include <cmath>
#include <map>
#include <tuple>

#include "uns/generate/HexMesh.h"

namespace uns
{
HexMesh make_cubed_sphere_shell(const CubedSphereSpec &s)
{
    const int M = s.ncube, nr = s.nr;
    HexMesh h;
    h.tag_names = {"", "inner", "outer"};
    // radial fractions along a ray: fluid layers (or all layers) grow geometrically
    const int nfluid = s.r_surf > 0 ? nr - s.n_solid : nr;
    std::vector<double> frac(nfluid + 1);
    {
        double sum = 0, dr = 1;
        std::vector<double> w(nfluid);
        for (int k = 0; k < nfluid; ++k) { w[k] = dr; sum += dr; dr *= s.r_stretch; }
        frac[0] = 0;
        for (int k = 0; k < nfluid; ++k) frac[k + 1] = frac[k] + w[k] / sum;
    }
    const Vec3 ec = 0.5 * (s.ell_lo + s.ell_hi), ea = 0.5 * (s.ell_hi - s.ell_lo);
    // outer radius along direction d : sphere, or intersection of the ray t*d with the ellipsoid
    auto outer_radius = [&](const Vec3 &d) {
        if (!s.ellipsoid) return s.r_out;
        // |((t d - c)/a)|^2 = 1
        const Vec3 dn(d.x / ea.x, d.y / ea.y, d.z / ea.z), cn(ec.x / ea.x, ec.y / ea.y, ec.z / ea.z);
        const double A = dot(dn, dn), Bq = -2 * dot(dn, cn), C = dot(cn, cn) - 1.0;
        return (-Bq + std::sqrt(Bq * Bq - 4 * A * C)) / (2 * A);
    };
    auto radius = [&](int k, const Vec3 &d) {
        if (s.r_surf > 0)
        {
            if (k <= s.n_solid) return s.r_in + (s.r_surf - s.r_in) * k / s.n_solid;
            return s.r_surf + (outer_radius(d) - s.r_surf) * frac[k - s.n_solid];
        }
        return s.r_in + (outer_radius(d) - s.r_in) * frac[k];
    };
    const double pi4 = std::atan(1.0);
    std::map<std::tuple<int, int, int, int>, Id> ids;
    // q in [0,M]^3 with one coordinate on the cube surface (q==0 or q==M)
    auto node = [&](const int q[3], int k) -> Id {
        const auto key = std::make_tuple(q[0], q[1], q[2], k);
        auto it = ids.find(key);
        if (it != ids.end()) return it->second;
        auto ang = [&](int i) { return std::tan(pi4 * (2.0 * i / M - 1.0)); }; // equiangular
        Vec3 d(ang(q[0]), ang(q[1]), ang(q[2]));
        d *= 1.0 / norm(d);
        const Id id = static_cast<Id>(h.x.size());
        h.x.push_back(radius(k, d) * d);
        ids.emplace(key, id);
        return id;
    };
    for (int axis = 0; axis < 3; ++axis)
        for (int side = 0; side <= 1; ++side)
            for (int iu = 0; iu < M; ++iu)
                for (int iv = 0; iv < M; ++iv)
                    for (int k = 0; k < nr; ++k)
                    {
                        std::array<Id, 8> v;
                        for (int q = 0; q < 8; ++q)
                        {
                            int lat[3];
                            lat[axis] = side ? M : 0;
                            lat[(axis + 1) % 3] = iu + (q & 1);
                            lat[(axis + 2) % 3] = iv + ((q >> 1) & 1);
                            v[q] = node(lat, k + ((q >> 2) & 1));
                        }
                        // local axes (u, v, radial) must be right-handed
                        const Vec3 e1 = h.x[v[1]] - h.x[v[0]], e2 = h.x[v[2]] - h.x[v[0]], e3 = h.x[v[4]] - h.x[v[0]];
                        if (dot(cross(e1, e2), e3) < 0) { std::swap(v[1], v[2]); std::swap(v[5], v[6]); }
                        h.hex.push_back(v);
                        h.material.push_back(k < s.n_solid ? kSolid : kFluid);
                        if (k == 0) h.tagged.push_back({{v[0], v[1], v[2], v[3]}, 1});
                        if (k == nr - 1) h.tagged.push_back({{v[4], v[5], v[6], v[7]}, 2});
                    }
    return h;
}
} // namespace uns
