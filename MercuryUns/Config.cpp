#include "Config.h"

#include <cmath>
#include <stdexcept>

namespace mercury
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
    c.lorentz_force = p.get_bool("lorentz_force", true);
    c.ideal.cuct = p.get_double("cuct", c.ideal.cuct);
    c.max_steps = static_cast<long>(p.get_double("max_Nstep", 1e6));
    c.t_end = p.get_double("max_Time", 1.0);
    c.log_every = p.get_int("log_every", 10);
    c.output_every = p.get_int("output_step", 0);
    c.output_dir = p.get_string("output_dir", "DATA");
    c.write_vtk = p.get_bool("write_vtk", true);
    c.write_fields = p.get_bool("write_fields", true);
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
        c.fsH = {rhoH, rhoH * Tn * c.sc.state_coeff_H};
        const double rhoNa = p.get_double("sw_rho_Na", 1e-6);
        c.fsNa = {rhoNa, rhoNa * Tn * c.sc.state_coeff_Na};
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

    c.bc_default = p.get_string("bc_default", "outflow");
    for (const std::string &k : p.keys())
        if (k.rfind("bc_", 0) == 0 && k != "bc_default") c.bc_by_tag[k.substr(3)] = p.get_string(k, "");
    if (!c.bc_by_tag.count("inner")) c.bc_by_tag["inner"] = "wall";
    if (!c.bc_by_tag.count("outer")) c.bc_by_tag["outer"] = "farfield";
    return c;
}
} // namespace mercury
