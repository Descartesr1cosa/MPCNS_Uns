#include "boundary/BoundaryInfo.h"

#include <stdexcept>

namespace mhd
{
using namespace uns;

BoundaryInfo BoundaryInfo::build(const LocalMesh &lm, const Geometry &geo, const Config &cfg)
{
    const CellComplex &cx = lm.mesh.cx;
    BoundaryInfo bi;
    bi.face_kind.assign(cx.n_faces(), FaceKind::Interior);
    bi.edge_mode.assign(cx.n_edges(), 0);
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const int nc = cx.f2c.count(f);
        if (nc == 2)
        {
            const bool s0 = cx.cell_material[cx.f2c.idx[cx.f2c.begin(f)]] == kSolid;
            const bool s1 = cx.cell_material[cx.f2c.idx[cx.f2c.begin(f) + 1]] == kSolid;
            bi.face_kind[f] = (s0 && s1) ? FaceKind::Skip : (s0 || s1) ? FaceKind::Wall : FaceKind::Interior;
            continue;
        }
        const Id c = cx.f2c.idx[cx.f2c.begin(f)];
        if (cx.cell_material[c] == kSolid) { bi.face_kind[f] = FaceKind::Skip; continue; }
        const int tag = cx.face_tag[f];
        std::string name = tag > 0 && tag < static_cast<int>(lm.mesh.tag_names.size()) ? lm.mesh.tag_names[tag] : "";
        auto it = cfg.bc_by_tag.find(name);
        const std::string type = it != cfg.bc_by_tag.end() ? it->second : cfg.bc_default;
        if (type == "wall") bi.face_kind[f] = FaceKind::Wall;
        else if (type == "inflow") bi.face_kind[f] = FaceKind::Inflow;
        else if (type == "outflow") bi.face_kind[f] = FaceKind::Outflow;
        else if (type == "farfield")
        {
            const Vec3 n_out = (cx.f2c.sgn[cx.f2c.begin(f)] / geo.face_area_mag[f]) * geo.face_area[f];
            bi.face_kind[f] = dot(cfg.u_inf, n_out) < 0 ? FaceKind::Inflow : FaceKind::Outflow;
        }
        else
            throw std::runtime_error("unknown boundary type '" + type + "' for face tag '" + name + "' (wall|farfield|inflow|outflow)");
    }
    for (Id f = 0; f < cx.n_faces(); ++f)
        if (bi.face_kind[f] == FaceKind::Inflow)
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) bi.edge_mode[cx.f2e.idx[k]] = 2;
    for (Id f = 0; f < cx.n_faces(); ++f)
        if (bi.face_kind[f] == FaceKind::Wall)
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) bi.edge_mode[cx.f2e.idx[k]] = 1;
    return bi;
}
} // namespace mhd
