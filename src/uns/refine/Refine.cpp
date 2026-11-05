#include "uns/refine/Refine.h"

#include <algorithm>
#include <stdexcept>

namespace uns
{
std::vector<int> compute_target_levels(const Mesh &coarse, const Geometry &g, const RefineSpec &spec, const CoarseSolution *sol)
{
    const Id nc = coarse.cx.n_cells();
    std::vector<int> lev(nc, 0);
    for (Id c = 0; c < nc; ++c)
    {
        lev[c] = spec.region_level(g.cell_center[c]);
        if (coarse.cx.cell_material[c] == kSolid) lev[c] = std::min(lev[c], spec.solid_max_level);
    }
    (void)sol; // flow-based indicators: milestone M3
    return lev;
}

RefineResult build_refined_mesh(const Mesh &, const std::vector<int> &)
{
    throw std::logic_error("refine: build_refined_mesh is milestone M3 (see PLAN.md)");
}
} // namespace uns
