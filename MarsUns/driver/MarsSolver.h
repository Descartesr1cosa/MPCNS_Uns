#pragma once
// Mars branch of MercuryUns: ten-ion Hall-MHD (Ar+, CO2+, CO+, C+, H+, N2+, N+, NO+, O+, O2+) on the signed cell complex.
// Ported from the structured-grid MarsZ5 solver (branch Mars of MPCNS_Mercury), see MarsUns/docs/Mars.md.
//
//   fluid      one five-equation Euler fluid per ion, Rusanov flux (first order as in MarsZ5, optional limited MUSCL)
//   sources    induce/gyration  rho (q B L / m_s U) (u_s - u+) x B,  w_s J x B,  - w_s grad p_e,  ion-neutral drag, heating/cooling
//   chemistry  cell-local coupled backward-Euler photochemical network (MarsChemistry), point-implicit momentum/energy losses
//   induction  Phi (induced face flux), B_add = IMF + crustal dipoles; E = ideal UCT + ambipolar (explicit) + Hall (semi-implicit,
//              the same GMRES machinery as the Mercury branch) ; B -= dt d1(E) so div B stays at round-off
//   boundaries inner "wall": photochemical-equilibrium Dirichlet ghost state (H+: 0.3 x solar-wind density), E = 0 on its edges;
//              outer "farfield": inflow = solar wind (H+) / floor state (the other ions), E = -u x B_imf; outflow = zero gradient
//   time       forward Euler (as MarsZ5), CFL = min(per-ion fluid CFL, total-plasma fast-magnetosonic CFL)
#include <mpi.h>

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "config/Config.h"
#include "boundary/BoundaryInfo.h"
#include "fluid/FluxRHS.h"
#include "fluid/Reconstruction.h"
#include "induction/EdgeEMF.h"
#include "induction/IdealEMF.h"
#include "induction/ImplicitInduction.h"
#include "physics/MarsAtmosphereModel.h"
#include "physics/MarsChemistry.h"
#include "physics/MarsDipoles.h"
#include "physics/MarsSpecies.h"
#include "uns/comm/Exchanger.h"
#include "uns/core/Params.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mars
{
using namespace mhd;
constexpr int NS = static_cast<int>(MARS::NumSpecies);

// Reference scales, switches and solar-wind state of the Mars case (nondimensional where noted).
struct MarsParams
{
    Config cfg;                      // generic configuration handed to the shared modules (boundary, ideal EMF, implicit induction)
    double kB = 1.380649e-23, q = 1.602176634e-19, mu0 = 1.256637061e-6, NA = 6.02214076e23;
    double L = 3.396e6, U = 3.27e4, T = 1.286e5, B = 3e-9, n = 4e6;          // reference scales (SI)
    double rho_ref = 0, p_ref = 0, time_ref = 0;
    double hall_scale = 0, pressure_scale = 0, mag_scale = 0;                 // nondimensional coefficients of E_hall, E_ambipolar, J x B
    double ne_floor = 1e-2;
    double rho_floor = 1e-10, p_floor = 1e-10;
    double sw_rho = 1, sw_p = 1;                                              // H+ solar wind
    uns::Vec3 sw_u, sw_B;
    double wall_T = 1500.0;
    bool hall = true, ambipolar = true, crustal = false;
    std::string dipole_file = "Dipoles";
    double cutoff_km = 0.0;
    MARS::AtmosphereConfig atm;
    static MarsParams from_params(const uns::Params &p);
};

class MarsSolver
{
public:
    MarsSolver(MPI_Comm comm, const uns::LocalMesh &lm, const MarsParams &prm);
    void initialize();   // collective: initial state (or init_fields), background field, atmosphere
    void run();          // time loop until t_end / max_steps with logging and output
    double step(double dt_max);

    struct Diag
    {
        double time = 0, dt = 0, div_b = 0;
        long step = 0;
        std::array<double, NS> mass{}, rho_min{}, rho_max{}, p_min{}, u_max{};
        double b_max = 0, j_max = 0;
        double dt_x[3] = {0, 0, 0};
        int dt_kind = -1;
    };
    Diag diagnostics();  // collective
    void write_vtk();
    void write_global_fields(const std::string &path);
    const std::vector<double> &U(int s) const { return U_[s]; }
    const std::vector<double> &phi() const { return Phi_; }
    const std::vector<double> &phi_add() const { return PhiAdd_; }
    double time() const { return time_; }
    long steps() const { return nstep_; }
    // test hooks (collective)
    void prepare();
    void apply_implicit_sources(double dt) { prepare(); implicit_sources_(dt); } // test hook: only the point-implicit source step
    std::vector<double> &U_mutable(int s) { return U_[s]; }
    double time_step_limit();
    const std::vector<double> &B_total_cell() const { return Btot_; }
    const std::vector<double> &electron_density() const { return ne_; }
    const std::vector<double> &ion_velocity() const { return uplus_; }

private:
    void build_background_();
    void build_atmosphere_();
    Prim wall_state_(uns::Id face, int s);
    double local_dt_() const;
    mutable struct DtLimit { double dt = 1e300; uns::Vec3 x; int kind = -1; } lim_; // cell and wave family that set the last time step (kind: ion index, or NS = fast mode)
    void fluid_and_sources_(double dt);
    void emf_(std::vector<double> &E);
    void implicit_sources_(double dt);
    void chemistry_(double dt);
    void floors_();
    void log_(const Diag &d) const;

    MPI_Comm comm_;
    int rank_ = 0;
    const uns::LocalMesh &lm_;
    MarsParams prm_;
    Config cfg_;
    uns::Geometry geo_;
    uns::comm::Exchanger ex_;
    EdgeEMF emf_asm_;
    BoundaryInfo bc_;
    Stencil stencil_;
    ImplicitInduction ind_;
    std::vector<char> active_;
    MARS::AtmosphereModel atmosphere_;
    Prim fs_[NS];
    GhostFn ghost_;
    std::unordered_map<uns::Id, std::array<Prim, NS>> wall_cache_;

    std::vector<double> U_[NS], W_[NS], grad_[NS], rhs_[NS];
    std::vector<double> Phi_, PhiAdd_, PhiTot_, BaddCell_;
    std::vector<double> Bind_, Btot_, J_, uplus_, ne_, ne_eff_, pe_, dpe_;
    std::vector<double> neutral_, photo_;   // NumNeutral / NumPhoto per local cell [cm^-3], [cm^-3 s^-1]
    // node -> fluid cells (ambipolar nodal pressure)
    std::vector<std::int64_t> nptr_;
    std::vector<uns::Id> ncell_;
    void ambipolar_emf_(double *E) const;
    double time_ = 0;
    long nstep_ = 0;
    long gmres_total_ = 0;
    bool initialized_ = false;
    std::vector<std::pair<double, std::string>> pvd_;
};
} // namespace mars
