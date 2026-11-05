#pragma once
// Agglomeration FAS multigrid for the fluid part of the pseudo-time (steady) iteration.  See docs/Multigrid.md for the mathematics.
//
//   level 0 : the fine cells (owned fluid cells)
//   level l : clusters of level-(l-1) clusters of the SAME rank, built by heavy-face matching (shared area vector as weight)
//   coarse faces : sum of the fine area vectors between two clusters; boundary faces : sum per (cluster, kind)
//   halo : clusters owned by other ranks that touch an owned cluster ("remote clusters").  They are exchanged through the fine-cell
//          Exchanger: every owned fine cell carries the value of its cluster, the ghost fine cells give the remote values back.
//
// The coarse operator is the first-order Rusanov flux on the clusters + the same source terms with frozen B, J, grad p_e.
// Only U (both species) is corrected; B keeps being advanced by the fine-level CT step (the coarse levels see it as frozen data).
#include <mpi.h>

#include <initializer_list>
#include <utility>
#include <vector>

#include "Config.h"
#include "boundary/BoundaryInfo.h"
#include "fluid/Reconstruction.h"
#include "uns/comm/Exchanger.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
struct MgLevel
{
    uns::Id n_own = 0, n_rem = 0;
    uns::Id n_tot() const { return n_own + n_rem; }
    std::vector<double> vol;                       // n_tot
    std::vector<std::int64_t> fptr;                // own cluster -> neighbours (own or remote), outward area vector
    std::vector<uns::Id> fnb;
    std::vector<uns::Vec3> fA;
    std::vector<std::int64_t> bptr;                // own cluster -> boundary faces grouped by kind
    std::vector<FaceKind> bkind;
    std::vector<uns::Vec3> bA;
    std::vector<uns::Id> cell_cluster;             // per fine local cell: cluster of this level, -1: none (solid / not covered)
    std::vector<uns::Id> rep_ghost;                // per remote cluster: one fine ghost cell of it
    std::vector<uns::Id> parent;                   // own cluster -> own cluster of the next level
    // frozen data (restricted from the fine level once per cycle), n_tot each
    std::vector<double> Btot, J, dpe, photo, vst, zero3;
    // states (5 per cluster, n_tot) and work arrays
    std::vector<double> U[2], U0[2], W[2], Ux[2];
    std::vector<double> g[2], rate[2], rr[2];      // forcing, rate, restricted rate (5 per own cluster)
    std::vector<double> dt;                        // own: step of the last pseudo step
};

class Multigrid
{
public:
    Multigrid(MPI_Comm comm, const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc, uns::comm::Exchanger &ex, const Config &cfg,
              const Prim fs[2]);
    int n_levels() const { return static_cast<int>(lv_.size()); }
    const MgLevel &level(int l) const { return lv_[l]; }
    // Global cluster count per level (collective, for the log).
    std::vector<long> global_sizes();

    // One FAS cycle.  Uf / ratef: fine conservative state and fine pseudo-time rate (U^{+}-U)/dt of both species (5 per local cell, valid on owned fluid cells);
    // frozen fine fields (3 per cell: Btot, J, dpe; 1 per cell: photo, vst).  Returns the prolonged coarse-grid correction in dU (5 per local cell, owned fluid cells).
    void cycle(const double *const Uf[2], const double *const ratef[2], const double *Btot, const double *J, const double *dpe, const double *photo,
               const double *vst, double *const dU[2]);

    // exposed for the tests
    void restrict_fine(int l, const double *const Uf[2], double *const U[2]) const;
    void flux_rhs(int l, const double *const W[2], double *const rhs[2]) const; // first-order flux part only, own clusters
    void exchange(int l, std::initializer_list<std::pair<double *, int>> arrays);
    double seconds_in_cycle() const { return t_cycle_; }

private:
    void build_level(int l, const std::vector<uns::Id> &own_id, uns::Id n_own);
    uns::Id match(int l, std::vector<uns::Id> &parent) const;
    void prim_all(int l, uns::Id n, const double *const U[2], double *const W[2]) const;
    void dt_level(int l);
    void source_rhs(int l, double *const rhs[2]);
    void pseudo_step(int l, bool use_forcing);                  // U^l <- Step(U^l) + dt g
    void rate_of(int l, bool with_forcing);                     // rate[l] = (Step(U^l) - U^l)/dt (+ g)
    void restrict_up(int l);                                    // level l -> l+1: U, rr
    void prolong_down(int l, double relax);                     // level l -> l-1 (l >= 2): U^{l-1} += relax * (U^l - U0^l)(parent)
    void floors(int l);

    MPI_Comm comm_;
    const uns::LocalMesh &lm_;
    const uns::Geometry &geo_;
    const BoundaryInfo &bc_;
    uns::comm::Exchanger &ex_;
    Config cfg_;
    Prim fs_[2];
    std::vector<MgLevel> lv_;
    std::vector<double> buf_;
    double t_cycle_ = 0;
    int nrank_ = 1;
};

// Largest alpha in {1, 1/2, 1/4, ...} (or 0) such that Uold + alpha dU keeps rho >= keep*rho_old and the internal energy >= keep*e_old
// (and, if grow > 0, rho <= grow*rho_old and e <= grow*e_old).
double positivity_safeguard(const double *Uold, const double *dU, double gamma, double keep = 0.1, double grow = 0.0);
} // namespace mercury
