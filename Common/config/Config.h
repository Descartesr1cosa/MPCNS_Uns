#pragma once
// Run configuration (nondimensional).  Built from a setup file by Config::from_params; tests fill it directly.
#include <map>
#include <string>

#include "induction/IdealEMF.h"
#include "physics/ChemistryRates.h"
#include "physics/Scales.h"
#include "uns/core/Params.h"
#include "uns/core/Types.h"

namespace mhd
{
struct SpeciesState
{
    double rho = 1, p = 1;
};

struct Config
{
    Scales sc;

    // ---- numerics
    double cfl = 1.0;                // dt = cfl V / sum_f (|u.n| + c_fast) A  (legacy definition)
    bool second_order = true;       // limited linear reconstruction (MUSCL); false = first order
    double recon_rel_limit = 0.0;   // theta of the relative (positivity) limiter on rho, p of the reconstruction, 0 = off (Stencil::set_rel_limit)
    double rho_floor = 1e-10, p_floor = 1e-12;
    double gyro_cfl = 0.0;          // optional dt <= gyro_cfl / Omega_max (not needed: the gyration is integrated exactly, physics/Gyration)
    double rho_floor_dt = 1e-2;     // density floor used in the wave speeds of the CFL estimate
    bool chemistry = false;         // Na+ photo-ionisation, drag, thermal relaxation (legacy source terms)
    ChemistryRates chem;
    bool ambipolar = false;         // electron-pressure force on the ions and ambipolar electric field
    int time_order = 2;             // 2: SSP-RK2 (explicit EMF, resistive split step); 1: forward Euler + implicit induction (E unknown)
    bool hall = false;              // Hall electric field (only with time_order = 1)
    double hall_r_min = 1.0, hall_r_max = 1.012; // radial taper of the Hall term (legacy r_min / r_max)
    double ind_rtol = 1e-4;         // implicit induction solve (GMRES)
    int ind_max_it = 100, ind_restart = 20;
    double ind_ilu_shift = 0.3;     // relative diagonal shift of the incomplete factorisation
    bool ind_ilu = false;           // block-Jacobi ILU(0) preconditioner (false: Jacobi)
    bool resistive = false;         // Mercury interior resistivity (implicit, window eta_max * shape(r))
    double eta_r_inner = 0.84, eta_r_outer = 1.04, eta_width = 0.02;
    double res_rtol = 1e-6;
    int res_max_it = 200;
    bool lorentz_force = true;      // J x B and induce source terms in the momentum/energy equations
    IdealParams ideal;              // UCT strength of the ideal EMF

    // ---- steady-state (pseudo-time) mode: local time stepping and optional multigrid (docs/Multigrid.md)
    bool steady = false;            // every cell advances with its own CFL step; B by B -= d1(dt_e E) with dt_e = min over the cells at the edge
    double steady_tol = 0.0;        // stop when the relative residual rates of H+ and B drop below this (0: run to max_Nstep)
    double lts_solid_dt = 1.0;      // pseudo step of edges without fluid cells (resistive solid), in units of the largest fluid step
    double lts_max_ratio = 0.0;     // cap of the cell step at ratio * (smallest step): 1 reproduces the global step (test), 0 = no cap
    bool implicit_ideal = false;    // steady mode: ideal induction (linear in B for frozen u) treated implicitly, uniform edge step lts_b_ratio * dt_min
    double lts_b_ratio = 0.0;       // cap of the EDGE step (induction) at ratio * (smallest step); 0 = none, 1 = B advances with the global step (diagnostic)
    double lts_smooth = 0.0;        // gradient limit of the pseudo step field: dt_c <= (1 + kappa) dt_neighbour (0: off).  Smooth dt keeps the CT term grad(dt) x E small
    bool debug_res = false;         // print the face with the largest relative B rate at every report
    int mg_levels = 0;              // 0: no multigrid; otherwise the number of levels incl. the fine one (agglomeration FAS, see multigrid/)
    int mg_pre = 1, mg_post = 1, mg_coarse = 2; // smoothing steps (fine level) before / after the coarse-grid correction, sweeps on the coarse levels
    int mg_cycle_every = 1;         // one coarse-grid correction every N fine pseudo steps
    long mg_start = 0;              // first pseudo step with a coarse-grid correction (the fine smoother alone carries the initial transient)
    double mg_growth = 2.0;         // a correction may raise rho and the internal energy of a cell by at most this factor (0: no cap)
    double mg_relax = 1.0;          // damping of the prolonged correction
    int mg_agglom_passes = 2;       // heavy-face matching passes per level (cluster size ~ 2^passes)
    int mg_min_clusters = 8;        // stop coarsening below this many clusters per rank

    // ---- run control
    long max_steps = 1000000;
    double t_end = 1.0;
    int log_every = 10;
    int output_every = 0;           // steps between ParaView/field dumps (0: only at the end)
    int restart_every = 0;          // steps between restart files (fields_final.unsb is rewritten), 0: only at the end
    long output_start = 0;          // first step eligible for periodic dumps (debugging)
    std::string output_dir = "DATA";
    bool write_vtk = true;
    bool write_profile = true;      // load_profile.json/.txt (per-rank timings and counts, fitted cost model) at the end of run()
    bool write_fields = true;       // global fields file (coarse solution / restart) at the end

    // ---- initial / free-stream state
    SpeciesState fsH, fsNa;         // free-stream density and pressure
    uns::Vec3 u_inf{-1, 0, 0};      // nondimensional velocity of the H+ free stream
    bool na_seed_at_rest = false;   // Na+ seed (and its far-field state) at rest instead of co-moving (legacy)
    uns::Vec3 u_Na() const { return na_seed_at_rest ? uns::Vec3() : u_inf; }
    uns::Vec3 B_imf{0, 0, 0};       // nondimensional IMF: part of the static background B_add (like the legacy code); B_induced starts at 0
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
} // namespace mhd
