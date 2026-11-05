#pragma once
// Run configuration (nondimensional).  Built from a setup file by Config::from_params; tests fill it directly.
#include <map>
#include <string>

#include "induction/IdealEMF.h"
#include "physics/Scales.h"
#include "uns/core/Params.h"
#include "uns/core/Types.h"

namespace mercury
{
struct SpeciesState
{
    double rho = 1, p = 1;
};

struct Config
{
    Scales sc;

    // ---- numerics
    double cfl = 0.6;
    bool second_order = true;       // limited linear reconstruction (MUSCL); false = first order
    double rho_floor = 1e-10, p_floor = 1e-12;
    double rho_floor_dt = 1e-2;     // density floor used in the wave speeds of the CFL estimate
    bool lorentz_force = true;      // J x B and induce source terms in the momentum/energy equations
    IdealParams ideal;              // UCT strength of the ideal EMF

    // ---- run control
    long max_steps = 1000000;
    double t_end = 1.0;
    int log_every = 10;
    int output_every = 0;           // steps between ParaView/field dumps (0: only at the end)
    std::string output_dir = "DATA";
    bool write_vtk = true;
    bool write_fields = true;       // global fields file (coarse solution / restart) at the end

    // ---- initial / free-stream state
    SpeciesState fsH, fsNa;         // free-stream density and pressure
    uns::Vec3 u_inf{-1, 0, 0};      // nondimensional velocity
    uns::Vec3 B_imf{0, 0, 0};       // nondimensional IMF (induced part, uniform)
    SpeciesState solidState{1e-6, 1e-8};
    std::string init_fields;        // global fields file (from mpcns_pre refine / a previous run), optional

    // ---- background dipole  A = m x (r - c) / (|r - c|^2 + soft^2)^{3/2}
    bool dipole = false;
    uns::Vec3 dipole_m{0, 0, 0}, dipole_c{0, 0, 0};
    double dipole_soft = 1e-3;

    // ---- boundary conditions per face-tag name: wall | farfield | inflow | outflow
    std::map<std::string, std::string> bc_by_tag;
    std::string bc_default = "outflow";

    static Config from_params(const uns::Params &p);
};
} // namespace mercury
