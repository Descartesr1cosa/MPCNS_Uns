#pragma once
// Non-conforming refinement of a hexahedral complex (see PLAN.md, section 4).
//
//   levels  : per coarse cell a minimum level (flow indicators) ...
//   regions : ... and manual regions evaluated per leaf inside the coarse cell
//   build_refined_mesh : octree per coarse hex; nodes are keyed by integer lattice points on the
//             canonical coarse entities, faces are the common refinement of the two sides, hanging
//             nodes become vertices of the coarse neighbours (polyhedra).
//   transfer_B : coarse face flux -> fine face flux with d2 B = 0 in every fine cell
//   prolong_cells : conservative, limited cell-field prolongation
#include <vector>

#include "uns/geometry/Geometry.h"
#include "uns/refine/RefineSpec.h"

namespace uns
{
struct CoarseSolution
{
    std::vector<double> rho, p;     // cell density / pressure (indicator inputs)
    std::vector<Vec3> B_cell, J_cell;
};

struct RefineOptions
{
    bool spherical_snap = false; // new nodes keep the radial interpolation of the coarse corners (shell geometry)
    int max_balance_iter = 32;
};

struct RefineResult
{
    Mesh fine;
    std::vector<int> cell_level;        // per fine cell
    std::vector<Id> parent_cell;        // coarse cell of every fine cell (fine cells are contiguous per coarse cell)
    std::vector<Id> face_coarse_face;   // coarse face a fine face lies on (-1: interior of a coarse cell)
    std::vector<Id> face_owner_cell;    // fine cell whose lattice frame defines the face (the smaller side)
    int balance_iterations = 0;
};

// Per-coarse-cell minimum level from the flow indicators (all zero if sol == nullptr); solid cap applied.
std::vector<int> compute_target_levels(const Mesh &coarse, const Geometry &g, const RefineSpec &spec, const CoarseSolution *sol);

RefineResult build_refined_mesh(const Mesh &coarse, const RefineSpec &spec, const std::vector<int> &cell_levels,
                                const RefineOptions &opt = RefineOptions());

struct TransferReport
{
    double max_rel_div = 0;   // max over fine cells of |d2 B| / sum|B_f|
    double max_tile_defect = 0;// max |sum of tile fluxes - coarse flux| / max|coarse flux|
};
// `coarse_flux`: per coarse face (coarse canonical orientation).  Returns per fine face (fine canonical orientation).
std::vector<double> transfer_B(const Mesh &coarse, const Geometry &cg, const std::vector<double> &coarse_flux,
                               const RefineResult &rr, const Geometry &fg, TransferReport *rep = nullptr);

// ncomp values per cell; limited gradient (Barth-Jespersen) with a volume-weighted conservative correction.
std::vector<double> prolong_cells(const Mesh &coarse, const Geometry &cg, const std::vector<double> &U, int ncomp,
                                  const RefineResult &rr, const Geometry &fg, bool limit = true);
} // namespace uns
