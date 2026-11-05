#pragma once
// The only place where the magnetic field changes:   B -= dt * d1(E).
// E is a 1-cochain on edges (line integrals).  Every physical EMF contribution
// (ideal, Hall, resistive, artificial) must enter through an edge array that was
// assembled by EdgeEMF + Exchanger::assemble_edges; B never sees anything else.
#include <vector>

#include "uns/parallel/LocalMesh.h"

namespace mercury
{
// Faces that are updated on this rank: all faces of owned cells.  All their edges
// are complete after edge assembly, so owner and alias copies of a face receive the
// identical update (same E, same canonical loop order).
std::vector<char> active_faces(const uns::CellComplex &cx);

void faraday_update(const uns::CellComplex &cx, const std::vector<char> &active, const double *E, double *B, double dt);
} // namespace mercury
