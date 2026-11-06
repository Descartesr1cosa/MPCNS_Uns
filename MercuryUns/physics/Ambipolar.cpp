#include "physics/Ambipolar.h"

#include <algorithm>

namespace mercury
{
using namespace mhd;
using namespace uns;

void electron_pressure_gradient(const LocalMesh &lm, const Geometry &geo, const double *WH, const double *WNa, double *dpe)
{
    const CellComplex &cx = lm.mesh.cx;
    std::fill(dpe, dpe + 3 * static_cast<std::size_t>(cx.n_cells()), 0.0);
    auto pe = [&](Id c) { return WH[5 * static_cast<std::size_t>(c) + 4] + WNa[5 * static_cast<std::size_t>(c) + 4]; };
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        Vec3 sum;
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            double pf = pe(c);
            if (cx.f2c.count(f) == 2)
            {
                const Id o = cx.f2c.idx[cx.f2c.begin(f)] == c ? cx.f2c.idx[cx.f2c.begin(f) + 1] : cx.f2c.idx[cx.f2c.begin(f)];
                if (cx.cell_material[o] == kFluid) pf = 0.5 * (pf + pe(o));
            }
            sum += (pf * static_cast<double>(cx.c2f.sgn[k])) * geo.face_area[f];
        }
        const double iv = 1.0 / geo.cell_volume[c];
        dpe[3 * c] = sum.x * iv; dpe[3 * c + 1] = sum.y * iv; dpe[3 * c + 2] = sum.z * iv;
    }
}

AmbipolarEMF::AmbipolarEMF(const LocalMesh &lm, const Geometry &geo) : lm_(lm), geo_(geo)
{
    const CellComplex &cx = lm.mesh.cx;
    nptr_.assign(cx.n_nodes() + 1, 0);
    std::vector<std::pair<Id, Id>> pr;
    std::vector<Id> cn;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        cx.cell_nodes(c, cn);
        for (Id v : cn) pr.emplace_back(v, c);
    }
    for (auto &p : pr) ++nptr_[p.first + 1];
    for (Id v = 0; v < cx.n_nodes(); ++v) nptr_[v + 1] += nptr_[v];
    ncell_.resize(pr.size());
    std::vector<std::int64_t> fill(nptr_.begin(), nptr_.end() - 1);
    for (auto &p : pr) ncell_[fill[p.first]++] = p.second;
}

void AmbipolarEMF::add(const Config &cfg, const double *WH, const double *WNa, double *E) const
{
    const CellComplex &cx = lm_.mesh.cx;
    std::vector<double> pn(cx.n_nodes(), 0.0);
    std::vector<char> have(cx.n_nodes(), 0);
    auto node_p = [&](Id v) {
        if (have[v]) return pn[v];
        have[v] = 1;
        double s = 0, w = 0;
        for (auto k = nptr_[v]; k < nptr_[v + 1]; ++k)
        {
            const Id c = ncell_[k];
            s += geo_.cell_volume[c] * (WH[5 * static_cast<std::size_t>(c) + 4] + WNa[5 * static_cast<std::size_t>(c) + 4]);
            w += geo_.cell_volume[c];
        }
        return pn[v] = w > 0 ? s / w : 0.0;
    };
    for (Id e = 0; e < cx.ent[1].n_owned; ++e)
    {
        // density averaged to the edge with the dual weights (fluid cells only)
        double rH = 0, rN = 0, wsum = 0;
        for (auto k = geo_.ec_ptr[e]; k < geo_.ec_ptr[e + 1]; ++k)
        {
            const Id c = geo_.ec_cell[k];
            if (cx.cell_material[c] == kSolid) continue;
            rH += geo_.ec_weight[k] * WH[5 * static_cast<std::size_t>(c)];
            rN += geo_.ec_weight[k] * WNa[5 * static_cast<std::size_t>(c)];
            wsum += geo_.ec_weight[k];
        }
        if (wsum <= 0) continue;
        const NumInfo ni = cfg.sc.numinfo(rH / wsum, rN / wsum);
        Id a, b;
        cx.edge_nodes(e, a, b);
        if (nptr_[a + 1] == nptr_[a] || nptr_[b + 1] == nptr_[b]) continue;
        E[e] -= cfg.sc.pressure_coef * (node_p(b) - node_p(a)) / ni.ne_eff;
    }
}
} // namespace mercury
