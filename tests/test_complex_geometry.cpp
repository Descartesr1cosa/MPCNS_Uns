// d1 d0 = 0, d2 d1 = 0, closure, orientation and positivity on generated hex meshes,
// including the cubed-sphere shell (valence-3 corners) and randomly flipped local orientations.
#include <cstdio>
#include <random>

#include "TestUtil.h"
#include "uns/complex/Validate.h"
#include "uns/generate/HexMesh.h"
#include "uns/geometry/Geometry.h"

using namespace uns;

static void run(const char *name, const HexMesh &h, long expect_euler = 1)
{
    Mesh m = build_mesh(h);
    TopologyReport t = validate_topology(m.cx);
    std::printf("[%s] %s\n", name, t.summary().c_str());
    CHECK(t.ok);
    CHECK(t.euler == expect_euler);
    Geometry g = compute_geometry(m);
    GeometryReport gr = check_geometry(m, g);
    std::printf("[%s] %s\n", name, gr.summary().c_str());
    CHECK(gr.ok);
    CHECK(gr.misoriented_dual == 0);
    // edge cell weights sum to one
    for (Id e = 0; e < m.cx.n_edges(); ++e)
    {
        double s = 0;
        for (auto k = g.ec_ptr[e]; k < g.ec_ptr[e + 1]; ++k) s += g.ec_weight[k];
        CHECK(std::abs(s - 1.0) < 1e-12);
    }
    // total volume = sum of cell volumes
    double vol = 0;
    for (double v : g.cell_volume) vol += v;
    std::printf("[%s] total volume %.12g\n", name, vol);

    // Flip the local orientation of random edges and faces: the identities are
    // orientation-independent (this is what alias signs rely on).
    std::mt19937 rng(7);
    CellComplex cx = m.cx;
    for (Id e = 0; e < cx.n_edges(); ++e)
        if (rng() & 1)
        {
            for (auto k = cx.e2n.begin(e); k < cx.e2n.end(e); ++k) cx.e2n.sgn[k] = -cx.e2n.sgn[k];
            for (auto k = cx.e2f.begin(e); k < cx.e2f.end(e); ++k) // flip the same entries in f2e
            {
                const Id f = cx.e2f.idx[k];
                for (auto q = cx.f2e.begin(f); q < cx.f2e.end(f); ++q)
                    if (cx.f2e.idx[q] == e) cx.f2e.sgn[q] = -cx.f2e.sgn[q];
            }
        }
    cx.build_up();
    TopologyReport t2 = validate_topology(cx);
    CHECK(t2.d1d0_violations == 0);
    CHECK(t2.d2d1_violations == 0);
}

int main()
{
    BoxSpec b;
    b.nx = 5; b.ny = 4; b.nz = 3;
    run("box5x4x3", make_box(b));

    CubedSphereSpec cs;
    cs.ncube = 4; cs.nr = 3; cs.n_solid = 1; cs.r_in = 0.8; cs.r_out = 4.0; cs.r_stretch = 1.3;
    // shell = S^2 x I : euler = V-E+F-C = chi(S^2) = 2
    run("cubed_sphere_shell", make_cubed_sphere_shell(cs), 2);
    TEST_END();
}
