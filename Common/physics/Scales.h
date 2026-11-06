#pragma once
// Nondimensional coefficients of the two-ion MHD model (same definitions as the legacy MercurySolver).
//
//   rho_s n_s = rho_s M_ref / M_s      nondimensional number density of species s
//   momentum source of species s:
//       induce_coeff * n_s (u_s - u+) x B   +   inver_MA2 * w_s (J x B),      w_s = n_s / n_e,eff
//   energy source = u_s . (momentum source),   u+ = chi_H u_H + chi_Na u_Na,  chi_s = n_s / n_e
//   Alfven speed^2 = inver_MA2 |B|^2 / rho
#include "uns/core/Params.h"

namespace mhd
{
struct NumInfo
{
    double nH = 0, nNa = 0, ne = 0, ne_eff = 0, chiH = 1, chiNa = 0, wH = 0, wNa = 0;
};

struct Scales
{
    double gamma = 5.0 / 3.0;
    double M_ref = 1, M_H = 1, M_Na = 23;       // molar masses (only ratios matter)
    double inver_MA2 = 1;                       // B_ref^2 / (mu0 rho_ref U_ref^2)
    double induce_coeff = 0;                    // q L_ref B_ref n_ref / (rho_ref U_ref)
    double state_coeff_H = 1, state_coeff_Na = 1; // p = rho T coeff   (T nondimensional)
    double ne_floor = 0.01;

    // chemistry / ambipolar / resistivity coefficients (legacy AddSourceToRHS_Fluid, AddAmbipolarEdgeEMF_, inver_Rem)
    double T_ref = 1, L_ref_m = 1;      // K, m
    double a1_Na = 0;                   // mass source per (cm^-3 s^-1) of Na+ production
    double a4 = 0;                      // drag: L_ref / U_ref
    double a5 = 0, a6 = 0;              // temperature relaxation / energy of the new ions
    double pressure_coef = 0;           // ambipolar: (rho_ref U_ref) / (q_e L_ref B_ref n_ref)
    double hall_coef = 0;               // B_ref M_ref / (U q L mu0 rho_ref N_A): E_Hall = hall_coef / n_e,eff * (J x B)
    double inver_Rem = 0;               // eta_max / (mu0 U_ref L_ref): nondimensional magnetic diffusivity

    // From the legacy-style physical reference values (U_ref, T_ref, n_ref, ... see MercuryUns/INPUT/setup.txt).
    static Scales from_params(const uns::Params &p);
    // Directly nondimensional (tests): unit masses, no induce term.
    static Scales simple(double gamma, double inver_MA2, double mass_ratio_Na = 1.0);
    NumInfo numinfo(double rhoH, double rhoNa) const;
};
} // namespace mhd
