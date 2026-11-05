#pragma once
#include <string>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns
{
struct TopologyReport
{
    bool ok = true;
    std::array<long, kNumDims> n{};
    long euler = 0;             // V - E + F - C
    long boundary_faces = 0;
    long d1d0_violations = 0;   // faces whose d1 d0 row is not zero
    long d2d1_violations = 0;   // cells whose d2 d1 row is not zero
    long bad_edges = 0, bad_loops = 0, bad_face_cells = 0;
    std::vector<std::string> messages;
    std::string summary() const;
};

// Pure combinatorial checks (no coordinates).  Valid for the global complex and
// for per-rank sub-complexes (an open ghost rim just shows up as boundary faces).
TopologyReport validate_topology(const CellComplex &cx);
} // namespace uns
