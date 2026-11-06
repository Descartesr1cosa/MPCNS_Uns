#pragma once
// Perfect-gas helpers for one species: conservative U = (rho, rho u, E), primitive W = (rho, u, v, w, p).
#include <algorithm>
#include <cmath>

namespace mhd
{
struct Prim
{
    double rho = 1, u = 0, v = 0, w = 0, p = 1;
};

inline Prim to_prim(const double *U, double gamma, double rho_floor, double p_floor)
{
    Prim W;
    W.rho = std::max(U[0], rho_floor);
    W.u = U[1] / W.rho; W.v = U[2] / W.rho; W.w = U[3] / W.rho;
    const double ke = 0.5 * W.rho * (W.u * W.u + W.v * W.v + W.w * W.w);
    W.p = std::max((gamma - 1.0) * (U[4] - ke), p_floor);
    return W;
}
inline void to_cons(const Prim &W, double gamma, double *U)
{
    U[0] = W.rho; U[1] = W.rho * W.u; U[2] = W.rho * W.v; U[3] = W.rho * W.w;
    U[4] = W.p / (gamma - 1.0) + 0.5 * W.rho * (W.u * W.u + W.v * W.v + W.w * W.w);
}
inline double sound(const Prim &W, double gamma) { return std::sqrt(gamma * W.p / W.rho); }

// Rusanov flux through a face of unit normal n (pointing from L to R); returns the flux per unit area.
inline void rusanov(const Prim &L, const Prim &R, const double n[3], double gamma, double *F)
{
    double UL[5], UR[5], FL[5], FR[5];
    to_cons(L, gamma, UL);
    to_cons(R, gamma, UR);
    auto phys = [&](const Prim &W, const double *U, double *Fx) {
        const double un = W.u * n[0] + W.v * n[1] + W.w * n[2];
        Fx[0] = W.rho * un;
        Fx[1] = W.rho * un * W.u + W.p * n[0];
        Fx[2] = W.rho * un * W.v + W.p * n[1];
        Fx[3] = W.rho * un * W.w + W.p * n[2];
        Fx[4] = un * (U[4] + W.p);
    };
    phys(L, UL, FL);
    phys(R, UR, FR);
    const double sL = std::abs(L.u * n[0] + L.v * n[1] + L.w * n[2]) + sound(L, gamma);
    const double sR = std::abs(R.u * n[0] + R.v * n[1] + R.w * n[2]) + sound(R, gamma);
    const double s = std::max(sL, sR);
    for (int m = 0; m < 5; ++m) F[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * s * (UR[m] - UL[m]);
}
} // namespace mhd
