#include "uns/refine/RefineSpec.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace uns
{
bool Region::contains(const Vec3 &p) const
{
    switch (kind)
    {
    case Kind::Box:
        return p.x >= a[0] && p.y >= a[1] && p.z >= a[2] && p.x <= a[3] && p.y <= a[4] && p.z <= a[5];
    case Kind::Sphere:
    {
        const Vec3 d = p - Vec3(a[0], a[1], a[2]);
        return norm(d) <= a[3];
    }
    case Kind::Shell:
    {
        const double r = norm(p);
        return r >= a[0] && r <= a[1];
    }
    case Kind::CylinderX: // axis parallel to x through (y0,z0)
    {
        const double dy = p.y - a[0], dz = p.z - a[1];
        return dy * dy + dz * dz <= a[2] * a[2] && p.x >= a[3] && p.x <= a[4];
    }
    }
    return false;
}

RefineSpec RefineSpec::from_params(const Params &p)
{
    RefineSpec s;
    s.max_level = p.get_int("refine_max_level", s.max_level);
    s.budget_cells = static_cast<long>(p.get_double("refine_budget_cells", 0));
    s.balance_2to1 = p.get_bool("refine_balance_2to1", true);
    s.solid_max_level = p.get_int("refine_solid_max_level", 0);
    for (const std::string &k : p.keys())
    {
        const bool reg = k.rfind("refine_region_", 0) == 0, ind = k.rfind("refine_indicator_", 0) == 0;
        if (!reg && !ind) continue;
        std::istringstream is(p.get_string(k, ""));
        std::string kind;
        is >> kind;
        std::vector<double> v;
        for (double x; is >> x;) v.push_back(x);
        if (reg)
        {
            Region r;
            std::size_t need = 0;
            if (kind == "box") { r.kind = Region::Kind::Box; need = 6; }
            else if (kind == "sphere") { r.kind = Region::Kind::Sphere; need = 4; }
            else if (kind == "shell") { r.kind = Region::Kind::Shell; need = 2; }
            else if (kind == "cylinder_x") { r.kind = Region::Kind::CylinderX; need = 5; }
            else throw std::runtime_error(k + ": unknown region kind '" + kind + "'");
            if (v.size() != need + 1) throw std::runtime_error(k + ": expected " + std::to_string(need) + " numbers and a level");
            r.level = static_cast<int>(v.back());
            v.pop_back();
            r.a = v;
            s.regions.push_back(r);
        }
        else
        {
            IndicatorRule r;
            if (kind == "density") r.kind = IndicatorRule::Kind::Density;
            else if (kind == "pressure") r.kind = IndicatorRule::Kind::Pressure;
            else if (kind == "current") r.kind = IndicatorRule::Kind::Current;
            else throw std::runtime_error(k + ": unknown indicator '" + kind + "'");
            if (v.size() != 2) throw std::runtime_error(k + ": expected 'threshold level'");
            r.threshold = v[0];
            r.level = static_cast<int>(v[1]);
            s.indicators.push_back(r);
        }
    }
    return s;
}

int RefineSpec::region_level(const Vec3 &p) const
{
    int lev = 0;
    for (const Region &r : regions)
        if (r.contains(p)) lev = std::max(lev, r.level);
    return std::min(lev, max_level);
}
} // namespace uns
