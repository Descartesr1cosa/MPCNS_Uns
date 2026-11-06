#pragma once
// Ideal-MHD edge EMF  E = -u x B  as a 1-cochain on edges, for edges of any valence
// (hexahedra, cubed-sphere corners, hanging-node polyhedra, MPI interfaces).
//
// Face candidates (Gardiner-Stone type, Cartesian form of the legacy Mercury scheme)
//     F_f   = 1/2 (F_L + F_R) - 1/2 a_f (B_R^ind - B_L^ind),   F_X = (u_X.n) B_X^tot - u_X B_f^tot,n
//     E_f   = - n x F_f                                         (tangential to the face)
// with n the canonical face normal pointing from cell L to cell R (L = cell whose c2f sign is +1) and
// a_f = max(|u_L.n|, |u_R.n|).  A face on the domain boundary uses its single cell (no dissipation).
//
// Edge value (line integral along the stored edge orientation):
//     E_e = sum_f w_{e,f} E_f . t_e                             central part, w from the dual face of e
//         + 2 C sum_f omega_f a_{e,f} (t_e . (m_f x n_f)) B_n,f   transverse upwind (UCT) part
// m_f is the unit in-face direction pointing away from the edge, a_{e,f} = max over cells at e of |u.m_f|,
// B_n,f = phi_ind,f / |A_f|.  omega_f are the weights w_{e,f} minimally corrected (per edge, 2x2 solve) so
// that the UCT part vanishes identically for a uniform B whatever the valence; on a regular hexahedral
// edge omega = w = 1/4 and the legacy formula is reproduced exactly (see tests/test_ideal_emf.cpp).
//
// Only the *owned* edges of the rank are evaluated (all their incident cells and faces are local);
// EdgeEMF::finish() then copies the value to the alias edges on other ranks.
#include "induction/EdgeEMF.h"

namespace mhd
{
struct IdealParams
{
    double cuct = 0.5; // strength of the transverse upwind part (legacy: 0.5)
};

struct IdealInputs
{
    const double *u_cell = nullptr;     // 3 per local cell (owned + ghost, current)
    const double *B_cell_tot = nullptr; // 3 per local cell
    const double *B_cell_ind = nullptr; // 3 per local cell
    const double *phi_tot = nullptr;    // per local face: total flux, canonical orientation
    const double *phi_ind = nullptr;    // per local face: evolved (induced) flux
};

// Adds the ideal EMF of every owned edge to E (call between EdgeEMF::begin and finish).
void add_ideal_emf(const uns::LocalMesh &lm, const uns::Geometry &geo, const IdealInputs &in, const IdealParams &prm, double *E);
} // namespace mhd
