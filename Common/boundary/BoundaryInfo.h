#pragma once
// Classification of faces / edges for boundary conditions, from face tags and cell materials.
//
//   Interior   fluid | fluid
//   Wall       fluid | solid interface, or a boundary face with bc "wall"   (slip wall; E = 0 on its edges)
//   Inflow     boundary face with bc "inflow", or "farfield" with u_inf . n_out < 0  (Dirichlet; E = 0 on its edges)
//   Outflow    boundary face with bc "outflow", or "farfield" with u_inf . n_out >= 0 (zero gradient; B free)
//   Skip       solid | solid, or a boundary face of a solid cell (nothing to solve)
//
// Edge modes: on the edges of Wall faces E = 0 (stationary perfect conductor, flux through the wall frozen); on the edges of
// Inflow faces E = -u_inf x B_imf . dl (convective field of the imposed free stream, so a uniform stream is preserved).
// Wall wins where both meet.
#include <vector>

#include "config/Config.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mhd
{
enum class FaceKind : char { Interior = 0, Wall, Inflow, Outflow, Skip };

struct BoundaryInfo
{
    std::vector<FaceKind> face_kind;  // per local face
    std::vector<char> edge_mode;      // per local edge: 0 free, 1 E = 0 (wall), 2 E = free-stream convective field (inflow)
    static BoundaryInfo build(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg);
    // outward-state of the ghost side of a boundary face for a given kind (see Euler helpers)
};
} // namespace mhd
