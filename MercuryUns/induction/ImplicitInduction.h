#pragma once
// Semi-implicit induction step with the edge EMF E as the unknown (first order in time).
//
//   Phi^{n+1} = Phi^n - dt d1(E)                                   -> div B is preserved for ANY E, i.e. for any solver accuracy
//   E         = E_expl + F(Phi^{n+1})                              E_expl: ideal + ambipolar EMF evaluated at t^n (and the boundary values)
//   F(Phi)    = M d1^T S2 Phi  +  H(Phi)                           implicit linear terms
//
//   resistive:  M = diag(eta_e / star1_e)  (Mercury interior window), d1^T S2 Phi = current flux through the dual face of e
//   Hall:       H(Phi)_e = sum_c w_{e,c} taper(r_c) alpha_c (J_c[Phi] x B0_c) . t_e ,   alpha_c = hall_coef / n_e,eff,c
//               B0 = total cell field at t^n is frozen (linearised Hall: B^{n+1} x J^{n+1} -> B^n x J^{n+1}); J_c[Phi] is the Gauss curl of the
//               cell field reconstructed from Phi.  This is the standard semi-implicit Hall treatment: unconditionally stable for whistlers.
//               The weights w_{e,c} are the dual-face shares of the cells at the edge, valid for any edge valence, and the edge sum goes through the
//               same owner reduction as every other EMF.
//   Substituting Phi^{n+1}:    (I + dt F d1) E = E_expl + F(Phi^n)
// The operator is not symmetric (Hall part ~ skew) with spectrum in the right half plane: restarted GMRES, Jacobi preconditioner (right).
// Dirichlet edges (imposed inflow EMF) are not unknowns.  Edges without implicit coupling are identity rows.
#include <vector>

#include "Config.h"
#include "induction/IdealEMF.h"
#include "uns/comm/Exchanger.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
class ImplicitInduction
{
public:
    ImplicitInduction(const uns::LocalMesh &lm, const uns::Geometry &geo, uns::comm::Exchanger &ex, const Config &cfg);
    bool active() const { return resistive_ || hall_ || ideal_; }
    bool has_ideal() const { return ideal_; }
    bool has_hall() const { return hall_; }

    // Frozen data of the step: total cell field B0 (3 per local cell) and electron density n_e,eff per local cell (0 for non-fluid cells).
    void set_frozen_state(const double *B0_cell, const double *ne_eff_cell);

    // Implicit ideal induction (cfg.implicit_ideal): the part of the ideal EMF that is linear in the induced flux for frozen u (central average of
    // u x B and the transverse upwind term, background field excluded) joins F.  `u_cell`: 3 per local cell (kept by pointer, must stay valid during solve),
    // `mask[e] != 0`: edges that carry the ideal EMF (free edges; walls / inflow edges excluded).
    void set_ideal_state(const double *u_cell, const std::vector<char> &mask);
    // L(phi') for face fluxes valid on all local faces; result valid on all local edges (assembled).
    void ideal_linear(const double *phi_all_faces, double *edge_out);

    // F(Phi') for face fluxes Phi' valid on ALL local faces; result valid on all local edges (assembled through the owners).
    void apply_F(const double *phi_all_faces, double *edge_out);

    // Solves for E (valid on all local edges).  `dirichlet[e] != 0`: E[e] = E_expl[e] is imposed.  Phi_n must be valid on all local faces.
    // Returns the number of GMRES iterations; `rel_res` receives the final relative residual.
    // `dt_edge` (optional, one value per local edge) replaces the scalar dt: local pseudo-time stepping, system (I + F d1 diag(dt_e)) E = ...
    int solve(const double *Phi_n, const double *E_expl, const std::vector<char> &dirichlet, double dt, double *E, double *rel_res = nullptr,
              const double *dt_edge = nullptr);

private:
    void apply_A(const std::vector<double> &x, std::vector<double> &y, double dt);
    double dot_owned(const std::vector<double> &a, const std::vector<double> &b);
    void build_precond_pattern();
    void factor_precond(double dt, const std::vector<char> &dir);
    void apply_precond(const std::vector<double> &v, std::vector<double> &z, const std::vector<char> &dir);

    const uns::LocalMesh &lm_;
    const uns::Geometry &geo_;
    uns::comm::Exchanger &ex_;
    bool resistive_, hall_, ideal_ = false;
    IdealParams ideal_prm_;
    const double *u_ideal_ = nullptr;
    std::vector<char> ideal_mask_;
    std::vector<double> ide_;
    double rtol_, hall_coef_ = 0;
    const double *dte_ = nullptr;             // per-edge pseudo-time step during solve() (null: scalar dt)
    int max_it_, restart_;

    std::vector<double> m_;                   // eta_e / star1_e (0 outside the resistive window)
    std::vector<double> taper_;               // Hall radial taper per local cell
    std::vector<double> alpha_;               // taper * hall_coef / ne_eff per cell (set_frozen_state)
    std::vector<double> B0_;                  // frozen total cell field
    std::vector<double> recon_;               // 3 coefficients per c2f entry: B_c = sum_k recon_k * Phi_{f_k}
    std::vector<std::int64_t> ce_ptr_;        // owned cell -> (edge, weight*t_e) list
    std::vector<uns::Id> ce_edge_;
    std::vector<uns::Vec3> ce_wt_;
    // block-Jacobi ILU(0) preconditioner of  I + dt M d1^T S2 d1  on the owned edges (couplings to other ranks dropped, Hall part ignored)
    bool use_ilu_ = true;
    double shift_ = 0.0;
    bool use_guess_ = true;
    std::vector<double> corr_;
    std::vector<std::int64_t> pc_ptr_;
    std::vector<uns::Id> pc_col_;
    std::vector<double> pc_coef_, pc_val_;
    std::vector<std::int64_t> pc_diag_;
    std::vector<double> diag_;                // Jacobi diagonal estimate
    std::vector<double> d2_;                  // sum_f star2 d1^2 over owned faces, assembled
    // work arrays
    std::vector<double> phi_, bcell_, jcell_, part_, tmp_;
};
} // namespace mercury
