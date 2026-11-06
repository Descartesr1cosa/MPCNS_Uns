#include "config/Config.h"

#include <cmath>
#include <stdexcept>

namespace mhd
{
using uns::Vec3;

Config Config::from_params(const uns::Params &p)
{
    Config c;
    c.sc = Scales::from_params(p);
    c.cfl = p.get_double("CFL", c.cfl);
    c.second_order = p.get_bool("second_order", true);
    c.rho_floor = p.get_double("rho_floor", c.rho_floor);
    c.p_floor = p.get_double("p_floor", c.p_floor);
    c.rho_floor_dt = p.get_double("rho_floor_dt", c.rho_floor_dt);
    c.gyro_cfl = p.get_double("gyro_cfl", c.gyro_cfl);
    c.lorentz_force = p.get_bool("lorentz_force", true);
    c.chemistry = p.get_bool("chemistry", false);
    c.chem.sk1 = p.get_double("sk1", c.chem.sk1);
    c.chem.sk2 = p.get_double("sk2", c.chem.sk2);
    c.chem.Tn0 = p.get_double("Tn0", c.chem.Tn0);
    c.ambipolar = p.get_bool("is_Ambipolar_Efield", false);
    c.resistive = p.get_bool("is_Mercury_resistance", false);
    c.time_order = p.get_int("time_order", 2);
    c.hall = p.get_bool("is_Hall_Efield", false);
    c.hall_r_min = p.get_double("r_min", c.hall_r_min);
    c.hall_r_max = p.get_double("r_max", c.hall_r_max);
    c.ind_rtol = p.get_double("induction_rtol", c.ind_rtol);
    c.ind_max_it = p.get_int("induction_max_it", c.ind_max_it);
    c.ind_restart = p.get_int("induction_restart", c.ind_restart);
    c.ind_ilu = p.get_bool("induction_ilu", false);
    c.ind_ilu_shift = p.get_double("induction_ilu_shift", c.ind_ilu_shift);
    c.eta_r_inner = p.get_double("mercury_eta_r_inner", c.eta_r_inner);
    c.eta_r_outer = p.get_double("mercury_eta_r_outer", c.eta_r_outer);
    c.eta_width = p.get_double("mercury_eta_width", c.eta_width);
    c.res_rtol = p.get_double("resistive_rtol", c.res_rtol);
    c.res_max_it = p.get_int("resistive_max_it", c.res_max_it);
    c.ideal.cuct = p.get_double("cuct", c.ideal.cuct);
    c.steady = p.get_bool("steady", false);
    c.steady_tol = p.get_double("steady_tol", 0.0);
    c.lts_solid_dt = p.get_double("lts_solid_dt", c.lts_solid_dt);
    c.lts_max_ratio = p.get_double("lts_max_ratio", 0.0);
    c.implicit_ideal = p.get_bool("implicit_ideal", false);
    c.lts_b_ratio = p.get_double("lts_b_ratio", 0.0);
    c.lts_smooth = p.get_double("lts_smooth", 0.0);
    c.debug_res = p.get_bool("debug_res", false);
    c.mg_levels = p.get_int("mg_levels", 0);
    c.mg_pre = p.get_int("mg_pre", c.mg_pre);
    c.mg_post = p.get_int("mg_post", c.mg_post);
    c.mg_coarse = p.get_int("mg_coarse", c.mg_coarse);
    c.mg_cycle_every = p.get_int("mg_cycle_every", c.mg_cycle_every);
    c.mg_start = static_cast<long>(p.get_double("mg_start", 0));
    c.mg_growth = p.get_double("mg_growth", c.mg_growth);
    c.mg_relax = p.get_double("mg_relax", c.mg_relax);
    c.mg_agglom_passes = p.get_int("mg_agglom_passes", c.mg_agglom_passes);
    c.mg_min_clusters = p.get_int("mg_min_clusters", c.mg_min_clusters);
    c.max_steps = static_cast<long>(p.get_double("max_Nstep", 1e6));
    c.t_end = p.get_double("max_Time", 1.0);
    c.log_every = p.get_int("log_every", 10);
    c.output_every = p.get_int("output_step", 0);
    c.restart_every = p.get_int("restart_step", 0);
    c.output_start = static_cast<long>(p.get_double("output_start", 0));
    c.output_dir = p.get_string("output_dir", "DATA");
    c.write_vtk = p.get_bool("write_vtk", true);
    c.write_fields = p.get_bool("write_fields", true);
    c.write_profile = p.get_bool("write_profile", true);
    c.init_fields = p.get_string("init_fields", "");

    // free stream: physical solar-wind values (like the legacy INITIAL list) or directly nondimensional
    if (p.has("sw_n"))
    {
        const double NA = p.get_double("NA", 6.02e23), U = p.get_double("U_ref", 4.031e5), T = p.get_double("T_ref", 1.7407e5);
        const double nref = p.get_double("n_ref", 4e7), B = p.get_double("B_ref", 1e-9);
        const double cy = p.get_double("sw_cy", 0.0), cz = p.get_double("sw_cz", 0.0);
        const double cx = -std::sqrt(std::max(0.0, 1.0 - cy * cy - cz * cz));
        const double Usw = p.get_double("sw_U", U);
        c.u_inf = {cx * Usw / U, cy * Usw / U, cz * Usw / U};
        const double rhoH = p.get_double("sw_n", nref) / nref * c.sc.M_H / c.sc.M_ref;
        const double Tn = p.get_double("sw_T", T) / T;
        const double rho_small = p.get_double("sw_rho_small", 1e-6); // Na+ seed (legacy): rho = rho_small rho_H, p = p_H rho_small / 23, at rest
        c.fsH = {rhoH, rhoH * Tn * c.sc.state_coeff_H};
        c.fsNa = {rho_small * rhoH, c.fsH.p * rho_small / 23.0};
        c.na_seed_at_rest = true;
        c.B_imf = {p.get_double("sw_Bx", 0) / B, p.get_double("sw_By", 0) / B, p.get_double("sw_Bz", 0) / B};
        (void)NA;
    }
    else
    {
        c.fsH = {p.get_double("init_rho_H", 1.0), p.get_double("init_p_H", 1.0)};
        c.fsNa = {p.get_double("init_rho_Na", 1e-6), p.get_double("init_p_Na", 1e-6)};
        const auto u = p.get_doubles("init_u"), b = p.get_doubles("init_B");
        if (u.size() == 3) c.u_inf = {u[0], u[1], u[2]};
        if (b.size() == 3) c.B_imf = {b[0], b[1], b[2]};
    }
    c.solidState = {p.get_double("solid_rho", c.solidState.rho), p.get_double("solid_p", c.solidState.p)};

    c.dipole = p.get_bool("dipole", false);
    const auto m = p.get_doubles("dipole_m"), cc = p.get_doubles("dipole_c");
    if (m.size() == 3) c.dipole_m = {m[0], m[1], m[2]};
    if (cc.size() == 3) c.dipole_c = {cc[0], cc[1], cc[2]};
    c.dipole_soft = p.get_double("dipole_soft", c.dipole_soft);
    // dipole in SI (legacy Dipole1: centre [m], moment [A m^2], softening [m]) -> nondimensional like the legacy DipoleField
    if (p.has("dipole_moment_SI"))
    {
        const auto mSI = p.get_doubles("dipole_moment_SI"), cSI = p.get_doubles("dipole_center_m");
        const double mu0 = p.get_double("mu_mag", 1.256637061e-6), B = p.get_double("B_ref", 1e-9), L = p.get_double("L_ref", 2.44e6);
        const double k = mu0 / (4.0 * 3.1415926535 * B * L * L * L);
        c.dipole = true;
        if (mSI.size() == 3) c.dipole_m = {mSI[0] * k, mSI[1] * k, mSI[2] * k};
        if (cSI.size() == 3) c.dipole_c = {cSI[0] / L, cSI[1] / L, cSI[2] / L};
        c.dipole_soft = p.get_double("dipole_a_m", 5.0e2) / L;
    }

    c.bc_default = p.get_string("bc_default", "outflow");
    for (const std::string &k : p.keys())
        if (k.rfind("bc_", 0) == 0 && k != "bc_default") c.bc_by_tag[k.substr(3)] = p.get_string(k, "");
    if (!c.bc_by_tag.count("inner")) c.bc_by_tag["inner"] = "wall";
    if (!c.bc_by_tag.count("outer")) c.bc_by_tag["outer"] = "farfield";
    return c;
}
} // namespace mhd
