#pragma once
// Hall electric field  E_hall = alpha (J x B)  on edges.   *** NOT IMPLEMENTED YET ***
//
// This is a deliberately empty hook (see docs/HallTodo.md for the agreed design):
//  - cell values  alpha_c J_c x B_c  are averaged onto each edge with the dual-face
//    shares w_{e,c} of uns::Geometry (valence independent) and enter EdgeEMF like the
//    ideal EMF, so they share the edge assembly and the owner;
//  - whistler-type diffusion with the same weights;
//  - h^2-tightened Hall time-step limit and sub-cycling / semi-implicit treatment.
// Until then the call is a no-op, so that the induction path (B -= dt d1 E) is complete
// and testable without Hall.
#include "induction/EdgeEMF.h"

namespace mercury
{
struct HallConfig
{
    bool enabled = false;
    double coef = 0.0; // hall_coef
};
void add_hall_emf(const HallConfig &cfg, EdgeEMF &emf, double *E);
// Hall limit dt <= C h^2 / (alpha |B|); returns +inf while unimplemented.
double hall_time_step_limit(const HallConfig &cfg);
} // namespace mercury
