#pragma once
// Metric quantities of a Mesh.  Kept apart from the topology: the complex never
// sees coordinates, everything metric (Hodge stars, area vectors, volumes,
// reconstruction weights) lives here.
//
// Conventions
//  * face centre  = area-weighted centroid of the triangle fan of the loop;
//    face area vector A_f follows the stored loop orientation (right-hand rule),
//    sum_f c2f[c,f] A_f = 0 for every closed cell.
//  * cell volume  = sum of tetrahedra (x_c, c_f, v_i, v_i+1) over the face-centre
//    triangulation of the cell surface (independent of x_c for a watertight surface).
//  * dual of an edge e inside cell c: the triangles (m_e, c_f, x_c) over the faces f
//    of c containing e.  Their vector areas are oriented along the edge direction.
//    w_{e,c} = |dual share of c| / |dual of e| (sums to 1 over the cells at e) is the
//    weight used to average cell quantities onto edges for any edge valence.
#include <string>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns
{
struct Geometry
{
    std::vector<Vec3> edge_vec, edge_mid;
    std::vector<double> edge_len;
    std::vector<Vec3> face_area, face_center;
    std::vector<double> face_area_mag;
    std::vector<Vec3> cell_center;
    std::vector<double> cell_volume, cell_min_tet_volume;

    std::vector<Vec3> edge_dual_area;   // vector area of the dual face of each edge
    std::vector<std::int64_t> ec_ptr;   // edge -> incident cells (CSR)
    std::vector<Id> ec_cell;
    std::vector<double> ec_weight;      // w_{e,c}
    std::vector<int> ec_misoriented;    // # dual pieces with negative projection on the edge (per edge)

    // Diagonal Hodge stars (all strictly positive for a valid mesh).
    std::vector<double> star0, star1, star2, star3;
};

struct GeometryReport
{
    bool ok = true;
    double min_volume = 0, max_volume = 0, min_tet = 0, min_edge = 0, max_edge = 0;
    double max_closure = 0;       // max_c |sum_f c2f A_f| / (sum_f |A_f|)
    long nonpositive_volume = 0, nonpositive_tet = 0, misoriented_dual = 0, nonpositive_star = 0;
    std::string summary() const;
};

Geometry compute_geometry(const Mesh &m);
GeometryReport check_geometry(const Mesh &m, const Geometry &g);
} // namespace uns
