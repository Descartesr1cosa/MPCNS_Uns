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
                    int species, const Prim &freestream, const double *W, const double *grad, double *rhs);

// Per-cell static data of the source terms (local numbering).
struct SourceFields
{
    const double *photo = nullptr; // Na+ production rate [cm^-3 s^-1]
    const double *vst = nullptr;   // drag / relaxation frequency [1/s] (photo-ionisation frequency of lit / shadowed region)
    const double *dpe = nullptr;   // grad p_e (3 per cell) or null
};

// Sources of species s (0 = H, 1 = Na) added to rhs on the owned fluid cells:
//   Lorentz:     inver_MA2 * w_s (J x B)  (the induce term is applied by physics/Gyration),  energy u_mid . (...), u_mid = u + dt S / (2 rho)
//                (the mid-point velocity makes the work equal the kinetic-energy change, so the force does not drain the thermal
//                energy in strongly magnetised, low-beta cells)
//   ambipolar:   - w_s grad p_e,  energy - w_s u_s . grad p_e
//   (dt_cell, if given, replaces dt per cell: pseudo-time stepping)
//   chemistry:   Na+ mass a1 q and internal energy a6 q Tn0 ; drag - a4 rho u vst (and kinetic energy) ;
//                temperature relaxation a5 n_s b1 and thermal drift a4 rho u^2 b2  (legacy AddSourceToRHS_Fluid)
void add_sources(const uns::LocalMesh &lm, const Config &cfg, int species, const double *WH, const double *WNa, const double *uplus,
                 const double *Btot, const double *J, const SourceFields &sf, double dt, double *rhs, const double *dt_cell = nullptr);

// The same on a plain list of n cells (multigrid levels): `material` (nullable) marks solid cells to skip, arrays indexed 0..n-1.
void add_sources_n(uns::Id n, const int *material, const Config &cfg, int species, const double *WH, const double *WNa, const double *uplus,
                   const double *Btot, const double *J, const SourceFields &sf, double dt, double *rhs, const double *dt_cell = nullptr);

// CFL time step of every owned fluid cell (0 for solid and ghost cells; the caller forwards it).  Used by the pseudo-time (steady) mode.
void cell_time_steps(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg, const double *WH, const double *WNa,
                     const double *Btot, double *dtc);

// CFL time step of the owned fluid cells (local value; caller reduces over ranks).
double local_time_step(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg, const double *WH, const double *WNa,
                       const double *Btot);
} // namespace mercury
