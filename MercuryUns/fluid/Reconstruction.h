#pragma once
// Second-order reconstruction on polyhedral cells: inverse-distance least-squares gradient of the primitive
// variables with Barth-Jespersen limiting.  Boundary faces contribute a virtual neighbour at the mirrored
// cell centre carrying the boundary ghost state, so cells at the boundary keep a full-rank stencil.
#include <vector>

#include "boundary/BoundaryInfo.h"
#include "physics/Euler.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
// Ghost-side primitive state of a boundary face (n_out = outward unit normal of the existing cell).
// Wall: H+ is reflected (slip wall), Na+ is absorbed (floor-state ghost: ions that reach the surface leave the fluid).
Prim ghost_state(FaceKind kind, int species, const Prim &inside, const uns::Vec3 &n_out, const Prim &freestream);

class Stencil
{
public:
    Stencil(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc);
    // limited gradients (5 primitives x 3) of the owned fluid cells; others are zero.  W: 5 per cell.
    void limited_gradients(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc, const double *W,
                           int species, const Prim &freestream, bool limit, double *grad) const;
private:
    std::vector<std::int64_t> ptr_;
    std::vector<uns::Id> nb_;          // neighbour cell, or -(face+1) for a virtual boundary neighbour
    std::vector<double> w_;
    std::vector<uns::Vec3> d_;
    std::vector<double> inv_;          // 9 per owned cell (symmetric inverse of sum w d d^T)
    std::vector<char> valid_;
};
} // namespace mercury
