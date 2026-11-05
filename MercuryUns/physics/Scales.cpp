#include "physics/Scales.h"

#include <algorithm>
#include <cmath>

namespace mercury
{
Scales Scales::from_params(const uns::Params &p)
{
    Scales s;
    s.gamma = p.get_double("gamma", 1.6667);
    const double NA = p.get_double("NA", 6.02e23), R = p.get_double("R_uni", 8.314);
    const double mu0 = p.get_double("mu_mag", 1.256637061e-6), qe = p.get_double("q_e", 1.602e-19);
    const double U = p.get_double("U_ref", 4.031e5), T = p.get_double("T_ref", 1.7407e5), n = p.get_double("n_ref", 4e7);
    const double B = p.get_double("B_ref", 1e-9), L = p.get_double("L_ref", 2.44e6);
    s.M_ref = p.get_double("Molecular_mass", 0.001);
    s.M_H = p.get_double("mole_mass1", 0.001);
    s.M_Na = p.get_double("mole_mass2", 0.023);
    const double rho_ref = s.M_ref * n / NA;
    s.inver_MA2 = B * B / (U * U * mu0 * rho_ref);
    s.induce_coeff = p.get_double("induce_scale", 1.0) * qe * L * B * n / (rho_ref * U); // induce_scale: experiments only
    s.state_coeff_H = R * T / (s.M_H * U * U);
    s.state_coeff_Na = R * T / (s.M_Na * U * U);
    s.ne_floor = p.get_double("ne_hall_floor", 0.01);
    const double kB = R / NA, m_Na = s.M_Na / NA;
    s.T_ref = T;
    s.L_ref_m = L;
    s.a1_Na = (L / U) * (m_Na / rho_ref) * 1e6;
    s.a4 = L / U;
    s.a5 = 3.0 * L * kB * n / (rho_ref * U * U * U);
    s.a6 = 1.0e6 * L * kB / (rho_ref * U * U * U * (s.gamma - 1.0));
    s.pressure_coef = rho_ref * U / (qe * L * B * n);
    s.hall_coef = B * s.M_ref / (U * qe * L * mu0 * rho_ref * NA);
    s.inver_Rem = p.get_double("eta_max_mercury", 0.0) / (mu0 * U * L);
    return s;
}

Scales Scales::simple(double gamma, double inver_MA2, double mass_ratio_Na)
{
    Scales s;
    s.gamma = gamma;
    s.M_ref = s.M_H = 1.0;
    s.M_Na = mass_ratio_Na;
    s.inver_MA2 = inver_MA2;
    s.induce_coeff = 0;
    s.state_coeff_H = s.state_coeff_Na = 1.0;
    return s;
}

NumInfo Scales::numinfo(double rhoH, double rhoNa) const
{
    NumInfo n;
    n.nH = std::max(rhoH, 0.0) / M_H * M_ref;
    n.nNa = std::max(rhoNa, 0.0) / M_Na * M_ref;
    n.ne = n.nH + n.nNa;
    n.ne_eff = std::sqrt(n.ne * n.ne + ne_floor * ne_floor);
    if (n.ne > 1e-300) { n.chiH = n.nH / n.ne; n.chiNa = n.nNa / n.ne; }
    else { n.chiH = 1.0; n.chiNa = 0.0; }
    n.wH = n.nH / n.ne_eff;
    n.wNa = n.nNa / n.ne_eff;
    return n;
}
} // namespace mercury
