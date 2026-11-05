#pragma once
// Exact integration of the "induce" force of the two-ion model,
//     du_s/dt = (induce n_s / rho_s) (u_s - u+) x B ,        u+ = chi_H u_H + chi_Na u_Na ,
// which is a linear rotation of the relative velocity w = u_H - u_Na about B:
//     dw/dt = Omega_rel w x b ,   Omega_rel = (a_H chi_Na + a_Na chi_H) |B| ,   a_s = induce M_ref / M_s
//     du_H = f_H dw ,  du_Na = -f_Na dw ,  f_H = a_H chi_Na / (a_H chi_Na + a_Na chi_H) ,  f_Na = 1 - f_H .
// Integrating it exactly (Strang-split around the explicit step) removes the gyro-frequency from the time-step limit
// (explicit Runge-Kutta amplifies a gyration with Omega dt ~ 0.2 slowly but surely: it blew up the polar cusp of the
// Mercury run once pick-up ions had loaded the field-line footprints).  The change of kinetic energy is added to the total
// energy, so the thermal energy is untouched.
#include "Config.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
// UH, UNa: 5 per cell; B: total cell field (3 per cell).  Acts on the owned fluid cells.
void gyrate(const uns::LocalMesh &lm, const Config &cfg, const double *Btot, double dt, double *UH, double *UNa);
} // namespace mercury
