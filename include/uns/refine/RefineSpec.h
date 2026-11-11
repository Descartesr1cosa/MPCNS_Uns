#pragma once
// What to refine: flow-based indicators, manually specified regions, or both.
// The target level of a cell is the maximum over all sources, clipped to max_level.
//
// setup keys (see case_template/Mercury/03_MercuryRefine/refine.txt):
//   refine_max_level N
//   refine_region_<i>  box x0 y0 z0 x1 y1 z1 level
//   refine_region_<i>  sphere cx cy cz r level
//   refine_region_<i>  shell r0 r1 level          (about the origin)
//   refine_region_<i>  cylinder_x y0 z0 r x0 x1 level
//   refine_indicator_<i> density|pressure|current  threshold level
//   refine_budget_cells N                          (optional global cell budget)
#include <string>
#include <vector>

#include "uns/core/Params.h"
#include "uns/core/Types.h"

namespace uns
{
struct Region
{
    enum class Kind { Box, Sphere, Shell, CylinderX } kind = Kind::Box;
    std::vector<double> a; // geometric parameters in the order of the setup line
    int level = 1;
    bool contains(const Vec3 &p) const;
};

struct IndicatorRule
{
    enum class Kind { Density, Pressure, Current } kind = Kind::Density;
    double threshold = 0.1; // h|grad rho|/rho, h|grad p|/p or h|J|/|B|
    int level = 1;
};

struct RefineSpec
{
    int max_level = 2;
    long budget_cells = 0; // 0 = no budget
    bool balance_2to1 = true;
    int solid_max_level = 0;
    std::vector<Region> regions;
    std::vector<IndicatorRule> indicators;

    static RefineSpec from_params(const Params &p);
    // level requested by the manual regions alone at point p
    int region_level(const Vec3 &p) const;
};
} // namespace uns
