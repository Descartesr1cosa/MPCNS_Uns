#include "physics/Gyration.h"

#include <cmath>

#include "physics/Euler.h"

namespace mercury
{
using namespace uns;

void gyrate(const LocalMesh &lm, const Config &cfg, const double *B, double dt, double *UH, double *UNa, const double *dt_cell)
{
    gyrate_n(lm.mesh.cx.ent[3].n_owned, lm.mesh.cx.cell_material.data(), cfg, B, dt, UH, UNa, dt_cell);
}

void gyrate_n(Id n, const int *material, const Config &cfg, const double *B, double dt, double *UH, double *UNa, const double *dt_cell)
{
    const Scales &sc = cfg.sc;
    if (sc.induce_coeff <= 0 || !cfg.lorentz_force) return;
    const double aH = sc.induce_coeff * sc.M_ref / sc.M_H, aN = sc.induce_coeff * sc.M_ref / sc.M_Na;
    for (Id c = 0; c < n; ++c)
    {
        if (material && material[c] == kSolid) continue;
        double *uh = UH + 5 * static_cast<std::size_t>(c), *un = UNa + 5 * static_cast<std::size_t>(c);
        const Prim wh = to_prim(uh, sc.gamma, cfg.rho_floor, cfg.p_floor), wn = to_prim(un, sc.gamma, cfg.rho_floor, cfg.p_floor);
        const NumInfo ni = sc.numinfo(uh[0], un[0]);
        const Vec3 b(B[3 * c], B[3 * c + 1], B[3 * c + 2]);
        const double bm = norm(b);
        const double denom = aH * ni.chiNa + aN * ni.chiH;
        if (bm <= 0 || denom <= 0) continue;
        const Vec3 bh = (1.0 / bm) * b;
        const Vec3 uH_(wh.u, wh.v, wh.w), uN_(wn.u, wn.v, wn.w);
        const Vec3 w = uH_ - uN_;
        const Vec3 w_par = dot(w, bh) * bh, w_perp = w - w_par;
        const double th = denom * bm * (dt_cell ? dt * dt_cell[c] : dt);
        const Vec3 dw = (std::cos(th) - 1.0) * w_perp + std::sin(th) * cross(w_perp, bh);
        const double fH = aH * ni.chiNa / denom, fN = aN * ni.chiH / denom;
        const Vec3 dH = fH * dw, dN = -fN * dw;
        // momentum and kinetic energy (thermal energy unchanged)
        const Vec3 uH2 = uH_ + dH, uN2 = uN_ + dN;
        auto set = [&](double *u, const Prim &w0, const Vec3 &v0, const Vec3 &v1) {
            const double rho = u[0];
            u[1] = rho * v1.x; u[2] = rho * v1.y; u[3] = rho * v1.z;
            u[4] += 0.5 * rho * (dot(v1, v1) - dot(v0, v0));
            (void)w0;
        };
        set(uh, wh, uH_, uH2);
        set(un, wn, uN_, uN2);
    }
}
} // namespace mercury
