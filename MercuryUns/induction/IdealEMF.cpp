#include "induction/IdealEMF.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mercury
{
using namespace uns;

namespace
{
inline Vec3 vec3(const double *a, Id i) { return {a[3 * static_cast<std::size_t>(i)], a[3 * static_cast<std::size_t>(i) + 1], a[3 * static_cast<std::size_t>(i) + 2]}; }

// Orthonormal basis (e1, e2) of the plane normal to the unit vector t.
void plane_basis(const Vec3 &t, Vec3 &e1, Vec3 &e2)
{
    const Vec3 a = std::abs(t.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    e1 = cross(t, a);
    e1 *= 1.0 / norm(e1);
    e2 = cross(t, e1);
}
} // namespace

void add_ideal_emf(const LocalMesh &lm, const Geometry &geo, const IdealInputs &in, const IdealParams &prm, double *E)
{
    const CellComplex &cx = lm.mesh.cx;

    // ---- face candidates, evaluated lazily for the faces around owned edges ----
    std::vector<Vec3> Ef(cx.n_faces());
    std::vector<char> have(cx.n_faces(), 0);
    auto face_candidate = [&](Id f) -> const Vec3 & {
        if (have[f]) return Ef[f];
        have[f] = 1;
        const double A = geo.face_area_mag[f];
        const Vec3 n = (1.0 / A) * geo.face_area[f];
        const double Bn_tot = in.phi_tot[f] / A;
        auto flux = [&](Id c) { // (u.n) B - u B_n
            const Vec3 u = vec3(in.u_cell, c);
            return dot(u, n) * vec3(in.B_cell_tot, c) - Bn_tot * u;
        };
        Vec3 F;
        const int nc = cx.f2c.count(f);
        if (nc == 2)
        {
            Id L = cx.f2c.idx[cx.f2c.begin(f)], R = cx.f2c.idx[cx.f2c.begin(f) + 1];
            if (cx.f2c.sgn[cx.f2c.begin(f)] < 0) std::swap(L, R); // n points from L (outward there) to R
            const double a = std::max(std::abs(dot(vec3(in.u_cell, L), n)), std::abs(dot(vec3(in.u_cell, R), n)));
            F = 0.5 * (flux(L) + flux(R)) - (0.5 * a) * (vec3(in.B_cell_ind, R) - vec3(in.B_cell_ind, L));
        }
        else
            F = flux(cx.f2c.idx[cx.f2c.begin(f)]);
        Ef[f] = -1.0 * cross(n, F);
        return Ef[f];
    };

    std::vector<double> w, aef;
    std::vector<Vec3> m, nrm;
    for (Id e = 0; e < cx.ent[1].n_owned; ++e)
    {
        const std::int64_t f0 = cx.e2f.begin(e), f1 = cx.e2f.end(e);
        const int nf = static_cast<int>(f1 - f0);
        if (nf == 0 || geo.edge_len[e] <= 0) continue;
        bool any_fluid = false; // edges of solid cells only carry no flow: E = 0
        for (std::int64_t k = geo.ec_ptr[e]; k < geo.ec_ptr[e + 1] && !any_fluid; ++k) any_fluid = cx.cell_material[geo.ec_cell[k]] == uns::kFluid;
        if (!any_fluid) continue;
        const Vec3 &t = geo.edge_vec[e];
        const double L = geo.edge_len[e];
        const Vec3 th = (1.0 / L) * t;

        // central part: weighted face candidates
        double Ec = 0;
        w.assign(nf, 0.0); m.assign(nf, Vec3()); nrm.assign(nf, Vec3()); aef.assign(nf, 0.0);
        for (int i = 0; i < nf; ++i)
        {
            const Id f = cx.e2f.idx[f0 + i];
            w[i] = geo.ef_weight[f0 + i];
            Ec += w[i] * dot(face_candidate(f), t);
            Vec3 v = geo.face_center[f] - geo.edge_mid[e];
            v -= dot(v, th) * th;
            const double nv = norm(v);
            m[i] = nv > 0 ? (1.0 / nv) * v : Vec3();
            nrm[i] = (1.0 / geo.face_area_mag[f]) * geo.face_area[f];
        }
        // transverse speeds: max over the cells at the edge of |u . m_f|
        for (std::int64_t k = geo.ec_ptr[e]; k < geo.ec_ptr[e + 1]; ++k)
        {
            const Vec3 u = vec3(in.u_cell, geo.ec_cell[k]);
            for (int i = 0; i < nf; ++i) aef[i] = std::max(aef[i], std::abs(dot(u, m[i])));
        }
        // Weights corrected so that a uniform B gives no upwind term:  sum_f omega_f q_f = 0, omega = w (1 - q.lambda).
        // The constraint lives in the 2-D plane normal to the edge (a component along the edge, present on warped faces,
        // would distort the weights by O(1)); solving a 2x2 system in an orthonormal basis of that plane keeps the
        // problem well conditioned, whereas a regularised 3x3 solve cancels catastrophically along the edge direction.
        Vec3 e1, e2;
        plane_basis(th, e1, e2);
        double M11 = 0, M12 = 0, M22 = 0, g1 = 0, g2 = 0;
        std::vector<double> qa(nf), qb(nf);
        for (int i = 0; i < nf; ++i)
        {
            const Vec3 qv = (aef[i] * dot(th, cross(m[i], nrm[i]))) * nrm[i];
            qa[i] = dot(qv, e1);
            qb[i] = dot(qv, e2);
            M11 += w[i] * qa[i] * qa[i];
            M12 += w[i] * qa[i] * qb[i];
            M22 += w[i] * qb[i] * qb[i];
            g1 += w[i] * qa[i];
            g2 += w[i] * qb[i];
        }
        const double tr = M11 + M22;
        const double eps = 1e-14 * tr; // only matters if all q are parallel (rank 1) or vanish
        const double det = (M11 + eps) * (M22 + eps) - M12 * M12;
        double l1 = 0, l2 = 0;
        if (tr > 0 && det > 0)
        {
            l1 = ((M22 + eps) * g1 - M12 * g2) / det;
            l2 = ((M11 + eps) * g2 - M12 * g1) / det;
        }
        double Eu = 0;
        for (int i = 0; i < nf; ++i)
        {
            const Id f = cx.e2f.idx[f0 + i];
            const double omega = w[i] * (1.0 - (qa[i] * l1 + qb[i] * l2));
            const double sig = dot(th, cross(m[i], nrm[i]));
            Eu += omega * aef[i] * sig * (in.phi_ind[f] / geo.face_area_mag[f]);
        }
        E[e] += Ec + 2.0 * prm.cuct * L * Eu;
    }
}
} // namespace mercury
