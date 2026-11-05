#include "uns/complex/Validate.h"

#include <sstream>

namespace uns
{
std::string TopologyReport::summary() const
{
    std::ostringstream s;
    s << "nodes=" << n[0] << " edges=" << n[1] << " faces=" << n[2] << " cells=" << n[3]
      << " euler(V-E+F-C)=" << euler << " boundary_faces=" << boundary_faces
      << " d1d0_bad=" << d1d0_violations << " d2d1_bad=" << d2d1_violations
      << " bad_edges=" << bad_edges << " bad_loops=" << bad_loops << " bad_face_cells=" << bad_face_cells
      << (ok ? "  [OK]" : "  [FAIL]");
    return s.str();
}

TopologyReport validate_topology(const CellComplex &cx)
{
    TopologyReport r;
    for (int d = 0; d < kNumDims; ++d) r.n[d] = cx.count(static_cast<Dim>(d));
    r.euler = r.n[0] - r.n[1] + r.n[2] - r.n[3];

    // every edge: exactly one start (-1) and one end (+1), distinct
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        int s = 0, t = 0;
        Id a = -1, b = -1;
        for (std::int64_t k = cx.e2n.begin(e); k < cx.e2n.end(e); ++k)
        {
            (cx.e2n.sgn[k] < 0 ? s : t)++;
            (cx.e2n.sgn[k] < 0 ? a : b) = cx.e2n.idx[k];
        }
        if (s != 1 || t != 1 || a == b || cx.e2n.count(e) != 2) ++r.bad_edges;
    }
    // every face loop closes: consecutive edges share the node where one ends and the next starts
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const int m = cx.f2e.count(f);
        bool ok = m >= 3;
        for (int i = 0; ok && i < m; ++i)
        {
            const std::int64_t k0 = cx.f2e.begin(f) + i, k1 = cx.f2e.begin(f) + (i + 1) % m;
            Id a0, b0, a1, b1;
            cx.edge_nodes(cx.f2e.idx[k0], a0, b0);
            cx.edge_nodes(cx.f2e.idx[k1], a1, b1);
            const Id exit0 = cx.f2e.sgn[k0] > 0 ? b0 : a0;
            const Id enter1 = cx.f2e.sgn[k1] > 0 ? a1 : b1;
            if (exit0 != enter1) ok = false;
        }
        if (!ok) ++r.bad_loops;
    }
    // faces: one boundary cell or two interior cells with opposite signs
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const int c = cx.f2c.count(f);
        if (c == 1) ++r.boundary_faces;
        else if (c == 2)
        {
            if (cx.f2c.sgn[cx.f2c.begin(f)] + cx.f2c.sgn[cx.f2c.begin(f) + 1] != 0) ++r.bad_face_cells;
        }
        else ++r.bad_face_cells;
    }
    r.d1d0_violations = count_nonzero_products(cx.f2e, cx.e2n);
    r.d2d1_violations = count_nonzero_products(cx.c2f, cx.f2e);
    r.ok = r.bad_edges == 0 && r.bad_loops == 0 && r.bad_face_cells == 0 && r.d1d0_violations == 0 &&
           r.d2d1_violations == 0;
    return r;
}
} // namespace uns
