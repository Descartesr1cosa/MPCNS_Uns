#pragma once
// Second-order reconstruction on polyhedral cells: inverse-distance least-squares gradient of the primitive
// variables with Barth-Jespersen limiting.  Boundary faces contribute a virtual neighbour at the mirrored
// cell centre carrying the boundary ghost state, so cells at the boundary keep a full-rank stencil.
#include <functional>
#include <vector>

#include "boundary/BoundaryInfo.h"
#include "physics/Euler.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mhd
{
// Ghost-side primitive state of a boundary face (n_out = outward unit normal of the existing cell).
// Wall: H+ is reflected (slip wall), Na+ is absorbed (floor-state ghost: ions that reach the surface leave the fluid).
Prim ghost_state(FaceKind kind, int species, const Prim &inside, const uns::Vec3 &n_out, const Prim &freestream);

// Optional override of the ghost state of a boundary face (used by the Mars branch for the photochemical-equilibrium wall):
// (face, species, inside state, outward unit normal) -> ghost state.  Null / empty = ghost_state() above.
using GhostFn = std::function<Prim(uns::Id face, int species, const Prim &inside, const uns::Vec3 &n_out)>;

class Stencil
{
public:
    Stencil(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc);
    // limited gradients (5 primitives x 3) of the owned fluid cells; others are zero.  W: 5 per cell.
    void limited_gradients(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc, const double *W,
                           int species, const Prim &freestream, bool limit, double *grad, const GhostFn *ghost = nullptr) const;
    // Positivity limiter (0 = off): after Barth-Jespersen the face values of rho and p may differ from the cell mean by at most
    // theta * mean, psi <- min(psi, theta |q_c| / |g.dx_f|).  Barth-Jespersen only bounds by the NEIGHBOUR extremes, so a trace cell next to a dense one
    // can be extrapolated tens of times above its own value and a forward-Euler step then makes rho or p negative (docs/MercuryRun.md, Diagnosis).
    void set_rel_limit(double theta) { rel_limit_ = theta; }
private:
    double rel_limit_ = 0.0;
    std::vector<std::int64_t> ptr_;
    std::vector<uns::Id> nb_;          // neighbour cell, or -(face+1) for a virtual boundary neighbour
    std::vector<double> w_;
    std::vector<uns::Vec3> d_;
    std::vector<double> inv_;          // 9 per owned cell (symmetric inverse of sum w d d^T)
    std::vector<char> valid_;
};
} // namespace mhd
