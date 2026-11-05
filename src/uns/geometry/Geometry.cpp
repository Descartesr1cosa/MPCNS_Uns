#include "uns/geometry/Geometry.h"

#include <algorithm>
#include <limits>
#include <sstream>

namespace uns
{
namespace
{
struct FaceTri
{
    Vec3 a, b; // v_i, v_{i+1}
};
} // namespace

Geometry compute_geometry(const Mesh &m)
{
    const CellComplex &cx = m.cx;
    const Id ne = cx.n_edges(), nf = cx.n_faces(), nc = cx.n_cells(), nn = cx.n_nodes();
    Geometry g;
    g.edge_vec.resize(ne); g.edge_mid.resize(ne); g.edge_len.resize(ne);
    for (Id e = 0; e < ne; ++e)
    {
        Id a, b;
        cx.edge_nodes(e, a, b);
        g.edge_vec[e] = m.x[b] - m.x[a];
        g.edge_mid[e] = 0.5 * (m.x[a] + m.x[b]);
        g.edge_len[e] = norm(g.edge_vec[e]);
    }

    g.face_area.resize(nf); g.face_center.resize(nf); g.face_area_mag.resize(nf);
    std::vector<Id> loop;
    for (Id f = 0; f < nf; ++f)
    {
        cx.face_loop(f, loop);
        const std::size_t k = loop.size();
        Vec3 p0;
        for (Id v : loop) p0 += m.x[v];
        p0 *= 1.0 / static_cast<double>(k);
        Vec3 cw;
        double wsum = 0;
        for (std::size_t i = 0; i < k; ++i)
        {
            const Vec3 &a = m.x[loop[i]], &b = m.x[loop[(i + 1) % k]];
            const double w = 0.5 * norm(cross(a - p0, b - p0));
            cw += w * ((p0 + a + b) * (1.0 / 3.0));
            wsum += w;
        }
        const Vec3 c = wsum > 0 ? cw * (1.0 / wsum) : p0;
        Vec3 A;
        for (std::size_t i = 0; i < k; ++i)
            A += 0.5 * cross(m.x[loop[i]] - c, m.x[loop[(i + 1) % k]] - c);
        g.face_center[f] = c;
        g.face_area[f] = A;
        g.face_area_mag[f] = norm(A);
    }

    g.cell_center.resize(nc); g.cell_volume.assign(nc, 0); g.cell_min_tet_volume.assign(nc, 0);
    g.star0.assign(nn, 0); g.star1.assign(ne, 0); g.star2.assign(nf, 0); g.star3.assign(nc, 0);
    g.edge_dual_area.assign(ne, Vec3());
    g.ec_misoriented.assign(ne, 0);

    // edge -> incident cells (sorted unique), weights filled below
    g.ec_ptr.assign(ne + 1, 0);
    std::vector<Id> tmp;
    for (Id e = 0; e < ne; ++e)
    {
        tmp.clear();
        for (std::int64_t k = cx.e2f.begin(e); k < cx.e2f.end(e); ++k)
        {
            const Id f = cx.e2f.idx[k];
            for (std::int64_t q = cx.f2c.begin(f); q < cx.f2c.end(f); ++q) tmp.push_back(cx.f2c.idx[q]);
        }
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        g.ec_cell.insert(g.ec_cell.end(), tmp.begin(), tmp.end());
        g.ec_ptr[e + 1] = static_cast<std::int64_t>(g.ec_cell.size());
    }
    g.ec_weight.assign(g.ec_cell.size(), 0.0);
    g.ef_weight.assign(cx.e2f.idx.size(), 0.0);
    std::vector<double> dual_len_face(nf, 0.0);

    std::vector<Id> cn;
    for (Id c = 0; c < nc; ++c)
    {
        cx.cell_nodes(c, cn);
        Vec3 xc0;
        for (Id v : cn) xc0 += m.x[v];
        xc0 *= 1.0 / static_cast<double>(cn.size());

        // pass 1: tet-volume weighted centroid with apex xc0
        double vsum = 0;
        Vec3 vc;
        for (std::int64_t kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            const Id f = cx.c2f.idx[kf];
            const int s = cx.c2f.sgn[kf];
            cx.face_loop(f, loop);
            const std::size_t k = loop.size();
            const Vec3 &cf = g.face_center[f];
            for (std::size_t i = 0; i < k; ++i)
            {
                const Vec3 &a = m.x[loop[i]], &b = m.x[loop[(i + 1) % k]];
                const Vec3 at = (0.5 * s) * cross(a - cf, b - cf);
                const double vt = dot(at, cf - xc0) / 3.0;
                vc += vt * ((xc0 + cf + a + b) * 0.25);
                vsum += vt;
            }
        }
        const Vec3 xc = vsum > 0 ? vc * (1.0 / vsum) : xc0;
        g.cell_center[c] = xc;

        // pass 2: volume, nodal dual volume, face-dual lengths, edge dual pieces
        double vol = 0, mintet = std::numeric_limits<double>::max();
        for (std::int64_t kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            const Id f = cx.c2f.idx[kf];
            const int s = cx.c2f.sgn[kf];
            cx.face_loop(f, loop);
            const std::size_t k = loop.size();
            const Vec3 &cf = g.face_center[f];
            dual_len_face[f] += norm(cf - xc);
            for (std::size_t i = 0; i < k; ++i)
            {
                const Id va = loop[i], vb = loop[(i + 1) % k];
                const Vec3 &a = m.x[va], &b = m.x[vb];
                const Vec3 at = (0.5 * s) * cross(a - cf, b - cf);
                const double vt = dot(at, cf - xc) / 3.0;
                vol += vt;
                mintet = std::min(mintet, vt);
                g.star0[va] += 0.5 * vt;
                g.star0[vb] += 0.5 * vt;

                // dual piece of edge (va,vb) in this cell via this face
                const std::int64_t ke = cx.f2e.begin(f) + static_cast<std::int64_t>(i);
                const Id e = cx.f2e.idx[ke];
                const int se = s * cx.f2e.sgn[ke]; // traversal along the stored edge direction?
                const Vec3 mid = g.edge_mid[e];
                const Vec3 piece = (0.5 * se) * cross(xc - mid, cf - mid);
                g.edge_dual_area[e] += piece;
                if (dot(piece, g.edge_vec[e]) <= 0) ++g.ec_misoriented[e];
                const auto lo = g.ec_cell.begin() + g.ec_ptr[e], hi = g.ec_cell.begin() + g.ec_ptr[e + 1];
                const auto it = std::lower_bound(lo, hi, c);
                g.ec_weight[it - g.ec_cell.begin()] += norm(piece);
                const auto flo = cx.e2f.idx.begin() + cx.e2f.begin(e), fhi = cx.e2f.idx.begin() + cx.e2f.end(e);
                g.ef_weight[std::lower_bound(flo, fhi, f) - cx.e2f.idx.begin()] += norm(piece);
            }
        }
        g.cell_volume[c] = vol;
        g.cell_min_tet_volume[c] = mintet;
        g.star3[c] = vol > 0 ? 1.0 / vol : 0.0;
    }

    for (Id e = 0; e < ne; ++e)
    {
        double s = 0;
        for (std::int64_t k = g.ec_ptr[e]; k < g.ec_ptr[e + 1]; ++k) s += g.ec_weight[k];
        if (s > 0)
        {
            for (std::int64_t k = g.ec_ptr[e]; k < g.ec_ptr[e + 1]; ++k) g.ec_weight[k] /= s;
            for (std::int64_t k = cx.e2f.begin(e); k < cx.e2f.end(e); ++k) g.ef_weight[k] /= s;
        }
        g.star1[e] = g.edge_len[e] > 0 ? norm(g.edge_dual_area[e]) / g.edge_len[e] : 0.0;
    }
    for (Id f = 0; f < nf; ++f)
        g.star2[f] = g.face_area_mag[f] > 0 ? dual_len_face[f] / g.face_area_mag[f] : 0.0;
    return g;
}

GeometryReport check_geometry(const Mesh &m, const Geometry &g)
{
    const CellComplex &cx = m.cx;
    GeometryReport r;
    r.min_volume = r.min_tet = r.min_edge = std::numeric_limits<double>::max();
    r.max_volume = r.max_edge = 0;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        Vec3 sum;
        double asum = 0;
        for (std::int64_t k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            sum += static_cast<double>(cx.c2f.sgn[k]) * g.face_area[cx.c2f.idx[k]];
            asum += g.face_area_mag[cx.c2f.idx[k]];
        }
        r.max_closure = std::max(r.max_closure, norm(sum) / std::max(asum, 1e-300));
        r.min_volume = std::min(r.min_volume, g.cell_volume[c]);
        r.max_volume = std::max(r.max_volume, g.cell_volume[c]);
        r.min_tet = std::min(r.min_tet, g.cell_min_tet_volume[c]);
        if (!(g.cell_volume[c] > 0)) ++r.nonpositive_volume;
        if (!(g.cell_min_tet_volume[c] > 0)) ++r.nonpositive_tet;
    }
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        r.min_edge = std::min(r.min_edge, g.edge_len[e]);
        r.max_edge = std::max(r.max_edge, g.edge_len[e]);
        if (g.ec_misoriented[e] > 0) ++r.misoriented_dual;
        if (!(g.star1[e] > 0)) ++r.nonpositive_star;
    }
    for (Id f = 0; f < cx.n_faces(); ++f) if (!(g.star2[f] > 0)) ++r.nonpositive_star;
    for (Id n = 0; n < cx.n_nodes(); ++n) if (!(g.star0[n] > 0)) ++r.nonpositive_star;
    for (Id c = 0; c < cx.n_cells(); ++c) if (!(g.star3[c] > 0)) ++r.nonpositive_star;
    r.ok = r.nonpositive_volume == 0 && r.nonpositive_star == 0 && r.max_closure < 1e-10;
    return r;
}

std::string GeometryReport::summary() const
{
    std::ostringstream s;
    s << "volume[min,max]=[" << min_volume << "," << max_volume << "] min_tet=" << min_tet
      << " edge[min,max]=[" << min_edge << "," << max_edge << "] max_closure=" << max_closure
      << " nonpos_vol=" << nonpositive_volume << " nonpos_tet=" << nonpositive_tet
      << " misoriented_dual_edges=" << misoriented_dual << " nonpos_star=" << nonpositive_star
      << (ok ? "  [OK]" : "  [FAIL]");
    return s.str();
}
} // namespace uns
