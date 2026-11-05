#pragma once
// Weighted cell-adjacency graph used for load balancing.
//   vertex = cell, weight = computational cost
//   edge   = shared face, weight = communication volume
#include <vector>

#include "uns/geometry/Geometry.h"

namespace uns
{
struct CostModel
{
    double fluid = 1.0;          // fluid + electromagnetics cell
    double solid = 0.35;         // electromagnetics-only (solid shell) cell
    double per_extra_face = 0.04;// relative cost of each face beyond 6 (polyhedra from hanging nodes)
    double comm_fluid = 1.0;     // halo volume of a face between two fluid cells (relative)
    double comm_solid = 0.3;     // halo volume of a face touching a solid cell
};

struct CellGraph
{
    std::vector<std::int64_t> xadj;
    std::vector<Id> adj;
    std::vector<double> adjw;
    std::vector<double> vw;
    std::vector<Vec3> centre; // for geometric partitioners
    Id n() const { return static_cast<Id>(vw.size()); }
};

CellGraph build_cell_graph(const Mesh &m, const Geometry &g, const CostModel &cost);
} // namespace uns
