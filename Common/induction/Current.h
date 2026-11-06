#pragma once
// Cell current  J = curl B  by Gauss' theorem on the polyhedron:
//     V_c J_c = sum_f (s_cf A_f) x B_f ,   B_f = 1/2 (B_L + B_R)  (one-sided on boundary faces).
// Only the tangential part of B_f enters, so the face-normal flux is not needed.
#include <vector>

#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mhd
{
// B: 3 per local cell; J (3 per local cell) is computed for the owned cells and zero elsewhere.
void cell_current(const uns::LocalMesh &lm, const uns::Geometry &geo, const double *B, double *J);
} // namespace mhd
