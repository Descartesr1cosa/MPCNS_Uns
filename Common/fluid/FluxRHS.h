#pragma once
// Finite-volume right-hand side of one ion fluid on polyhedral cells (model independent):
//   dU_c/dt = - 1/V_c sum_f s_cf F_f |A_f|     (Rusanov flux of linearly reconstructed face states)
#include "boundary/BoundaryInfo.h"
#include "fluid/Reconstruction.h"

namespace mhd
{
// rhs (5 per cell) is overwritten on the owned fluid cells and zero elsewhere.
void fluid_flux_rhs(const uns::LocalMesh &lm, const uns::Geometry &geo, const BoundaryInfo &bc, const Config &cfg,
                    int species, const Prim &freestream, const double *W, const double *grad, double *rhs,
                    const GhostFn *ghost = nullptr);
} // namespace mhd
