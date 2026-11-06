#pragma once
// Mercury interior resistivity, backward Euler, as a constrained-transport update.
//
//   dB/dt = - curl(eta J),   J = curl B (nondimensional, mu0 absorbed),   eta(r) = inver_Rem * shape(r),
//   shape = 1/2 [tanh((r - r_inner)/w) - tanh((r - r_outer)/w)]   (window r_inner..r_outer, smooth, legacy values 0.84, 1.04, 0.02)
//
// DEC form.  With D = d1 (faces x edges), S2 = diag(star2) and the edge EMF E_res = M D^T S2 Phi, M = diag(eta_e / star1_e):
//       Phi^{n+1} = Phi^n - dt D E_res ,   (M^{-1} + dt D^T S2 D) E_res = D^T S2 Phi^n                 (SPD in the edge unknowns)
// so B still changes only through d1 of an edge 1-cochain and div B is untouched.  Edges with eta = 0 carry E_res = 0.
// The system is solved by Jacobi-preconditioned CG in parallel: every face is processed by its owner rank, edge sums go through
// the owner (reduce + broadcast), dot products over owned edges only.  This also gives the E = eta J tangential EMF on the
// planet-surface wall edges (their ideal EMF is zero).
#include <mpi.h>

#include <vector>

#include "config/Config.h"
#include "uns/comm/Exchanger.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
using namespace mhd;
class ResistiveSolver
{
public:
    ResistiveSolver(const uns::LocalMesh &lm, const uns::Geometry &geo, uns::comm::Exchanger &ex, const Config &cfg);
    bool enabled() const { return enabled_; }
    long n_owned_edges() const { return n_owned_res_; }   // size of the local unknown set (load estimate)
    // Phi: face fluxes (valid on the owned faces); updated on the faces of owned cells.  Returns the CG iteration count.
    int advance(double *Phi, double dt, const std::vector<char> &active_faces);
    int last_iterations() const { return last_it_; }
private:
    void apply(const std::vector<double> &x, std::vector<double> &y, double dt);   // y = (M^-1 + dt D^T S2 D) x on R
    void dt_s2_d(const double *edge_vec, double *out) const;                     // out_e = sum_{owned f} d1 star2 (D x)_f  (assembled)
    double dot(const std::vector<double> &a, const std::vector<double> &b);

    const uns::LocalMesh &lm_;
    const uns::Geometry &geo_;
    uns::comm::Exchanger &ex_;
    double rtol_;
    int max_it_;
    bool enabled_ = false;
    std::vector<char> in_r_;
    std::vector<double> minv_;     // star1 / eta on R
    std::vector<double> d2_;       // sum_f star2 d1^2 (diagonal of D^T S2 D), assembled
    std::vector<double> x_, b_, r_, z_, p_, q_;
    long n_owned_res_ = 0;
    int last_it_ = 0;
};
} // namespace mercury
