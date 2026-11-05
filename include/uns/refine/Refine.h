#pragma once
// Refinement pipeline (milestone M3).  Interfaces are fixed here so that Preprocess
// (explicit `task=refine` step) and the solver start-up (library call) share them.
//
//   1. compute_target_levels : indicators (from a coarse solution) + manual regions
//                              -> level per coarse cell (<= max_level, 2:1 balanced)
//   2. build_refined_mesh    : lattice-key based non-conforming refinement; hanging nodes
//                              become vertices of the coarse neighbours (polyhedra)
//   3. transfer_B            : coarse face flux -> fine faces with  d2 B = 0  in every cell
//   4. transfer_U            : conservative cell prolongation
#include <vector>

#include "uns/geometry/Geometry.h"
#include "uns/refine/RefineSpec.h"

namespace uns
{
struct CoarseSolution
{
    std::vector<double> rho, p;   // cell density / pressure (indicator inputs)
    std::vector<Vec3> B_cell, J_cell;
};

struct RefineResult
{
    Mesh fine;
    std::vector<int> cell_level;           // level of every fine cell
    std::vector<Id> parent_cell;           // coarse cell of every fine cell
    std::vector<Id> face_parent;           // coarse face a fine face lies on (-1: interior of a coarse cell)
};

std::vector<int> compute_target_levels(const Mesh &coarse, const Geometry &g, const RefineSpec &spec,
                                       const CoarseSolution *sol /* may be null: regions only */);
// Not implemented yet: throws std::logic_error("refine: milestone M3").
RefineResult build_refined_mesh(const Mesh &coarse, const std::vector<int> &levels);
} // namespace uns
