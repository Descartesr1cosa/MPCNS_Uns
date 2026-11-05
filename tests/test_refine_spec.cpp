// Manual refinement regions and their parsing (flow-based indicators arrive with M3).
#include <cstdio>

#include "TestUtil.h"
#include "uns/generate/HexMesh.h"
#include "uns/refine/Refine.h"

using namespace uns;

int main()
{
    Params p;
    p.set_override("refine_max_level=3");
    p.set_override("refine_region_1=box -6 -2 -2 -1 2 2 3");
    p.set_override("refine_region_2=sphere 0 0 0 1.5 1");
    p.set_override("refine_region_3=shell 4 5 2");
    p.set_override("refine_indicator_1=current 0.3 2");
    RefineSpec s = RefineSpec::from_params(p);
    CHECK(s.regions.size() == 3 && s.indicators.size() == 1 && s.max_level == 3);
    CHECK(s.region_level({-3, 0, 0}) == 3);
    CHECK(s.region_level({0, 0, 1}) == 1);
    CHECK(s.region_level({4.5, 0, 0}) == 2);
    CHECK(s.region_level({10, 10, 10}) == 0);

    CubedSphereSpec cs;
    cs.ncube = 6; cs.nr = 8; cs.n_solid = 2; cs.r_in = 1; cs.r_out = 10; cs.r_stretch = 1.2;
    Mesh m = build_mesh(make_cubed_sphere_shell(cs));
    Geometry g = compute_geometry(m);
    s.solid_max_level = 0;
    auto lev = compute_target_levels(m, g, s, nullptr);
    long refined = 0;
    for (Id c = 0; c < m.cx.n_cells(); ++c)
    {
        refined += lev[c] > 0;
        if (m.cx.cell_material[c] == kSolid) CHECK(lev[c] == 0);
    }
    std::printf("cells=%d flagged=%ld\n", m.cx.n_cells(), refined);
    CHECK(refined == 0); // regions act at leaf level inside build_refined_mesh; indicators need a solution
    TEST_END();
}
