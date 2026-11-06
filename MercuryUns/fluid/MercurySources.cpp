#include "fluid/MercurySources.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mercury
{
using namespace mhd;
using namespace uns;

void add_sources(const LocalMesh &lm, const Config &cfg, int s, const double *WH, const double *WNa, const double *up, const double *B,
                 const double *J, const SourceFields &sf, double dt, double *rhs, const double *dt_cell)
{
    const CellComplex &cx = lm.mesh.cx;
    add_sources_n(cx.ent[3].n_owned, cx.cell_material.data(), cfg, s, WH, WNa, up, B, J, sf, dt, rhs, dt_cell);
}

void add_sources_n(Id n, const int *material, const Config &cfg, int s, const double *WH, const double *WNa, const double *up, const double *B,
                   const double *J, const SourceFields &sf, double dt, double *rhs, const double *dt_cell)
{
    const Scales &sc = cfg.sc;
    const double *W = s == 0 ? WH : WNa;
    for (Id c = 0; c < n; ++c)
    {
        if (material && material[c] == kSolid) continue;
        const std::size_t i = c;
        const double dtc = dt_cell ? dt_cell[i] : dt;
        const NumInfo ni = sc.numinfo(WH[5 * i], WNa[5 * i]);
        const double ns = s == 0 ? ni.nH : ni.nNa, ws = s == 0 ? ni.wH : ni.wNa;
        const double rho = std::max(W[5 * i], 0.0);
        const Vec3 u(W[5 * i + 1], W[5 * i + 2], W[5 * i + 3]);
        double *r = rhs + 5 * i;
        Vec3 S;
        if (cfg.lorentz_force)
        {
            const Vec3 uplus(up[3 * i], up[3 * i + 1], up[3 * i + 2]), b(B[3 * i], B[3 * i + 1], B[3 * i + 2]), j(J[3 * i], J[3 * i + 1], J[3 * i + 2]);
            (void)uplus; (void)ns; // the induce term is integrated exactly by physics/Gyration (operator splitting)
            S += (sc.inver_MA2 * ws) * cross(j, b);
        }
        if (cfg.ambipolar && sf.dpe) S -= ws * Vec3(sf.dpe[3 * i], sf.dpe[3 * i + 1], sf.dpe[3 * i + 2]);
        r[1] += S.x; r[2] += S.y; r[3] += S.z;
        r[4] += dot(u, S) + (rho > 1e-30 ? 0.5 * dtc * dot(S, S) / rho : 0.0); // work with the mid-point velocity
        if (cfg.chemistry && sf.vst)
        {
            const double vst = sf.vst[i], u2 = dot(u, u);
            const double Ts = W[5 * i + 4] / (rho > 0 ? rho : 1e-300) / (s == 0 ? sc.state_coeff_H : sc.state_coeff_Na) * sc.T_ref;
            const double Msum = sc.M_H + sc.M_Na;
            const double b1 = vst * (cfg.chem.Tn0 - Ts) * (s == 0 ? sc.M_H / Msum : 0.5);
            const double b2 = (s == 0 ? sc.M_Na / Msum : 0.5) * vst;
            r[1] -= sc.a4 * rho * u.x * vst; r[2] -= sc.a4 * rho * u.y * vst; r[3] -= sc.a4 * rho * u.z * vst;
            r[4] -= sc.a4 * rho * u2 * vst;               // drag
            r[4] += sc.a5 * ns * b1;                      // temperature relaxation
            r[4] += sc.a4 * rho * u2 * b2;                // thermal drift
            if (s == 1 && sf.photo)
            {
                r[0] += sc.a1_Na * sf.photo[i];           // Na+ creation
                r[4] += sc.a6 * sf.photo[i] * cfg.chem.Tn0; // internal energy of the new ions
            }
        }
    }
}

void cell_time_steps(const LocalMesh &lm, const Geometry &geo, const Config &cfg, const double *WH, const double *WNa, const double *B, double *dtc)
{
    const CellComplex &cx = lm.mesh.cx;
    const double inf = std::numeric_limits<double>::infinity();
    for (Id c = 0; c < cx.n_cells(); ++c) dtc[c] = 0.0;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        double dt = inf;
        const Vec3 b(B[3 * c], B[3 * c + 1], B[3 * c + 2]);
        const double B2 = dot(b, b);
        // dt <= cfl V / sum_f (|u_s.n| + c_s,fast) A for each species.  The Lorentz force acts on the ions in proportion to their
        // share w_s, so the magnetic wave speed is collective: it is evaluated with the total mass density (floored as in the legacy
        // code), while each species keeps its own sound speed.
        // explicit gyration of the ions about B: dt <= gyro_cfl / Omega (the lighter species gyrates fastest)
        if (cfg.gyro_cfl > 0 && cfg.sc.induce_coeff > 0)
        {
            const double omega = cfg.sc.induce_coeff * std::sqrt(B2) * std::max(cfg.sc.M_ref / cfg.sc.M_H, cfg.sc.M_ref / cfg.sc.M_Na);
            if (omega > 0) dt = std::min(dt, cfg.gyro_cfl / omega);
        }
        const double rho_tot = std::max(WH[5 * static_cast<std::size_t>(c)] + WNa[5 * static_cast<std::size_t>(c)], cfg.rho_floor_dt);
        for (int s = 0; s < 2; ++s)
        {
            const double *w = (s == 0 ? WH : WNa) + 5 * static_cast<std::size_t>(c);
            const double cs2 = cfg.sc.gamma * w[4] / std::max(w[0], 1e-30);
            double denom = 0;
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
            {
                const Id f = cx.c2f.idx[k];
                const double A = geo.face_area_mag[f];
                const Vec3 n = (1.0 / A) * geo.face_area[f];
                const double vA2 = cfg.sc.inver_MA2 * B2 / rho_tot;
                const double vAn2 = cfg.sc.inver_MA2 * dot(b, n) * dot(b, n) / rho_tot;
                const double term = cs2 + vA2;
                const double cf = std::sqrt(0.5 * (term + std::sqrt(std::max(0.0, term * term - 4 * cs2 * vAn2))));
                denom += (std::abs(w[1] * n.x + w[2] * n.y + w[3] * n.z) + cf) * A;
            }
            if (denom > 0) dt = std::min(dt, cfg.cfl * geo.cell_volume[c] / denom);
        }
        dtc[c] = dt;
    }
}

double local_time_step(const LocalMesh &lm, const Geometry &geo, const Config &cfg, const double *WH, const double *WNa, const double *B)
{
    std::vector<double> d(lm.mesh.cx.n_cells());
    cell_time_steps(lm, geo, cfg, WH, WNa, B, d.data());
    double dt = std::numeric_limits<double>::infinity();
    for (Id c = 0; c < lm.mesh.cx.ent[3].n_owned; ++c)
        if (lm.mesh.cx.cell_material[c] != kSolid) dt = std::min(dt, d[c]);
    return dt;
}
} // namespace mercury
