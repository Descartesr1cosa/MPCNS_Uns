#pragma once
// Two-ion MHD solver on the signed cell complex: second-order finite volume for the fluid, constrained transport
// for B (B -= dt d1(E), E an edge 1-cochain assembled through the owner), SSP-RK2 in time.
//
// One stage:  forward halo (U cells, B faces) -> derived fields (primitives, cell B, u+, J) -> limited gradients
//             (forwarded) -> fluid fluxes + sources -> ideal EMF on owned edges -> EdgeEMF::finish (copy to aliases).
// The Hall term is not part of this version (hall/HallEMF.h is a no-op hook).
#include <mpi.h>

#include <memory>
#include <string>
#include <vector>

#include "Config.h"
#include "boundary/BoundaryInfo.h"
#include "fluid/FluidRHS.h"
#include "fluid/Reconstruction.h"
#include "hall/HallEMF.h"
#include "induction/EdgeEMF.h"
#include "induction/IdealEMF.h"
#include "induction/ImplicitInduction.h"
#include "induction/Resistive.h"
#include "multigrid/Multigrid.h"
#include "physics/Ambipolar.h"
#include "uns/comm/Exchanger.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
struct State
{
    std::vector<double> U[2]; // H+, Na+ conservative variables, 5 per local cell
    std::vector<double> Phi;  // induced magnetic flux per local face (canonical orientation)
};

struct Diagnostics
{
    double time = 0, dt = 0;
    long step = 0;
    double div_b = 0;          // max over cells of |d2 B| / sum|B_f|  (global)
    double min_rho = 0, min_p = 0; // global minima over fluid cells (H+)
    double mass_H = 0, mass_Na = 0;
    double res_H = 0, res_Na = 0, res_B = 0, res_Bs = 0; // relative residual rates of the last step: |dU[0]/dt| (L2, volume weighted) / |U[0]|, and |dPhi/dt| / |Phi_tot|
    double ind_it = 0, ind_res = 0;          // GMRES iterations per step since the last report, worst relative residual since then
    double dt_min = 0, dt_max = 0;           // pseudo-time mode: smallest / largest cell step of the last step
    double dt_hall = 0;        // explicit whistler limit  min_c h_c^2 n_e / (pi^2 hall_coef |B|)  over the Hall-active cells (information only: the Hall term is implicit)
};

struct CostModelRecommendation
{
    double cost_fluid = 1, cost_solid = 0.35, comm_fluid = 1, comm_solid = 0.3;
};

class Solver
{
public:
    Solver(MPI_Comm comm, const uns::LocalMesh &lm, const Config &cfg);
    void initialize();              // initial state / background / boundary data (collective)
    void run();                     // time loop until t_end or max_steps, with logging and output
    double step(double dt_max);     // dispatches on Config::steady / Config::time_order
    double step_pseudo();           // local pseudo-time step (steady mode): cell steps dtc, edge steps dte = min over the cells at the edge
    double step_first_order(double dt_max);  // forward Euler fluid + semi-implicit induction (Hall, resistivity) with E as the unknown     // one SSP-RK2 step with the CFL step (<= dt_max); returns the dt used
    Diagnostics diagnostics();      // collective
    void write_vtk();               // collective: one .vtu per rank, .pvtu and .pvd
    // Collective: gathers the per-rank timings / entity counts of the time loop, fits a fluid/solid/communication cost model by
    // least squares and writes <base>.json (full record) and <base>.txt (cost_* keys readable by `mpcns_pre cost_profile=`).
    void write_load_profile(const std::string &base);
    void write_global_fields(const std::string &path); // collective: coarse solution / restart file (rank 0 writes)

    // wall-clock bookkeeping of the time loop (for load-balance calibration)
    struct Profile
    {
        double wall = 0, comm = 0;             // seconds spent in step(): total / in exchanges, reductions and waits
        long steps = 0, cg_iterations = 0;
    long n_ghost_cells = 0;
        long n_fluid = 0, n_solid = 0, n_edges = 0, n_res_edges = 0, n_faces = 0;   // owned entities of this rank
        double compute() const { return wall - comm; }
    };
    Profile profile() const;
    void reset_profile();

    const State &state() const { return s_; }
    State &state() { return s_; }
    const uns::LocalMesh &mesh() const { return lm_; }
    const uns::Geometry &geometry() const { return geo_; }
    double time() const { return time_; }
    long steps() const { return nstep_; }
    // derived cell fields of the last prepared state (3 per cell / 5 per cell)
    const std::vector<double> &prim(int species) const { return W_[species]; }
    const std::vector<double> &B_induced_cell() const { return Bind_; }
    const std::vector<double> &B_total_cell() const { return Btot_; }
    const std::vector<double> &current_cell() const { return J_; }

private:
    struct Rhs { std::vector<double> dU[2]; std::vector<double> E; };
    void prepare(State &s);                      // halo exchange + derived fields
    double local_dt() const;
    void rhs(State &s, Rhs &r, bool with_emf = true); // requires prepare(s); with_emf=false: fluid + sources only
    void advance(const State &base, const Rhs &r, double dt, State &out);
    void enforce_floors(State &s) const;
    Diagnostics collect(double dt);
    void log(const Diagnostics &d) const;

    MPI_Comm comm_;
    int rank_ = 0;
    const uns::LocalMesh &lm_;
    Config cfg_;
    uns::Geometry geo_;
    uns::comm::Exchanger ex_;
    EdgeEMF emf_;
    AmbipolarEMF amb_;
    ResistiveSolver res_;
    ImplicitInduction ind_;
    BoundaryInfo bc_;
    Stencil stencil_;
    std::vector<char> active_;
    Prim fs_[2];

    State s_, s1_, s2_;
    Rhs r0_, r1_;
    std::vector<double> PhiAdd_, BaddCell_;
    std::vector<double> W_[2], grad_[2], uplus_, Bind_, Btot_, J_, PhiTot_, dpe_, photo_, vst_;
    double wall_s_ = 0, dt_cur_ = 0;
    long cg_total_ = 0, cg_last_ = 0, step_last_ = 0;
    double ind_res_max_ = 0;
    std::vector<double> dtc_, dte_;           // pseudo-time: per-cell / per-edge step
    double dtc_min_ = 0, dtc_max_ = 0;
    double res_loc_[8] = {0, 0, 0, 0, 0, 0, 0, 0};  // local sums of the last step: sum V (dU0_H/dt)^2, sum V U0_H^2, same for Na, sum (dPhi/dt)^2, sum Phi^2
    void mg_correct_();                         // one FAS coarse-grid correction of U (steady mode, mg_levels > 1)
    std::unique_ptr<Multigrid> mg_;
    double mg_wall_ = 0;
    void smooth_dt_(int sweeps);                // gradient limiter of dtc_ (lts_smooth), forwards dtc_
    std::vector<double> dts_;                   // smoothed step field of the previous iteration (persistent)
    struct Worst { uns::Id f = -1; double rel = 0, circ = 0, phi = 0, dt = 0, pind = 0, dphi = 0; uns::Vec3 x; } worst_;
    void accumulate_B_residual_(const std::vector<double> &before, const double *dt_edge, double dt_uniform);
    void accumulate_residual_(const double *E); // after rhs(s_, r0_): residual rates of the step at U^n (before the update)
public:
    double multigrid_seconds() const { return mg_wall_; }
    const Multigrid *multigrid() const { return mg_.get(); }
    double last_induction_residual() const { return ind_res_max_; }
private:
    double time_ = 0;
    long nstep_ = 0;
    std::vector<std::pair<double, std::string>> pvd_;
    bool initialized_ = false;
};
} // namespace mercury
