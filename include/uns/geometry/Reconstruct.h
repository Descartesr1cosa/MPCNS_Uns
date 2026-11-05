#pragma once
// Cell-centred reconstruction of face fluxes and inverse-distance least-squares gradients.
#include <vector>

#include "uns/geometry/Geometry.h"

namespace uns
{
// Area-weighted least-squares vector b of every cell with  b . n_f = flux_f / A_f  (outward n_f);
// exact for uniform fields on any polyhedron with three independent face normals, and the
// average of opposite faces on parallelepipeds.  `flux` is in the canonical face orientation.
std::vector<Vec3> reconstruct_cell_vector(const Mesh &m, const Geometry &g, const std::vector<double> &flux);
// Same for the cells [0, ncells) with a raw pointer (local arrays).
std::vector<Vec3> reconstruct_cell_vector(const CellComplex &cx, const Geometry &g, const double *flux, Id ncells);
} // namespace uns
