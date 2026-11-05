#pragma once
// Finite-volume right-hand side of one species on polyhedral cells:
//   dU_c/dt = - 1/V_c sum_f s_cf F_f |A_f|     (Rusanov flux of linearly reconstructed face states)
// plus the electromagnetic momentum/energy sources of the two-ion model.
#include "boundary/BoundaryInfo.h"
#include "fluid/Reconstruction.h"

namespace mercury
{
// rhs (5 per cell) is overwritten on the owned fluid cells and zero elsewhere.
void fluid_flux_rhs(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc, const Config &cfg,
                    const Prim &freestream, const double *W, const double *grad, double *rhs);

// Lorentz-type sources for species s (0 = H, 1 = Na):
//   S_mom = induce * n_s (u_s - u+) x B + inver_MA2 * w_s (J x B),  S_E = u_s . S_mom      (added to rhs)
void add_em_sources(const uns::LocalMesh &lm, const Config &cfg, int species, const double *WH, const double *WNa, const double *uplus,
                    const double *Btot, const double *J, double *rhs);

// CFL time step of the owned fluid cells (local value; caller reduces over ranks).
double local_time_step(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg, const double *WH, const double *WNa,
                       const double *Btot);
} // namespace mercury
