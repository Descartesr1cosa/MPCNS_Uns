# MercuryUns — two-ion MHD solver (simple second-order version)

Code: `MercuryUns/` (driver `driver/Solver.cpp`).  Setup: `MercuryUns/INPUT/setup.txt`.  Tests: `tests/test_solver.cpp`, `tests/test_ideal_emf.cpp`.

## Model (nondimensional, same scaling as the legacy MercurySolver)
Species s = H+, Na+ with U_s = (rho, rho u, E), p = (gamma-1)(E - rho u²/2):

    dU_s/dt + div F(U_s) = S_s
    S_s,mom = induce * n_s (u_s - u+) x B  +  inver_MA2 * w_s (J x B),      S_s,E = u_s . S_s,mom
    n_s = rho_s M_ref / M_s ,  chi_s = n_s / n_e ,  w_s = n_s / n_e,eff ,  u+ = chi_H u_H + chi_Na u_Na
    dB/dt = - curl E ,   B = B_induced (evolved) + B_dipole (static)

Not included yet: photo-ionisation, drag/charge exchange, ambipolar electric field, resistivity, **Hall EMF** (hook only), the
thermal relaxation terms of the legacy source routine.

## Numerics
- **Fluid**: cell-centred finite volume on polyhedral cells; Rusanov flux; second order by linear reconstruction of the primitive variables
  (inverse-distance least-squares gradient over the face neighbours + Barth–Jespersen limiter, `fluid/Reconstruction.*`).  Boundary and
  fluid/solid faces contribute a mirrored virtual neighbour so boundary cells keep a full-rank stencil.  `second_order=0` switches to first order.
- **Magnetic field**: constrained transport only: `Phi -= dt d1(E)`; `E` = ideal EMF (`docs/IdealEMF.md`) assembled on the owned edges and copied to
  the alias edges.  The dipole is the exactly divergence-free flux `d1(A.dl)` of its vector potential.  `J = curl B_induced` per cell by Gauss' theorem.
- **Time**: SSP-RK2 for U and Phi alike (`Phi^{n+1} = Phi^n - dt/2 d1(E^n + E^{(1)})`, so div B stays at round-off).  `dt = CFL * 2V / sum_f (|u.n| + c_fast) A`
  per cell and species (density floored at `rho_floor_dt = 1e-2` in the wave speeds as in the legacy code), global minimum.
- **Parallel**: each stage forwards U (cells), Phi (faces) and the limited gradients (cells) from the owners, then evaluates; the EMF of an owned edge
  uses only local data (vertex-neighbour ghost layer).  Results agree with the serial run to round-off (see tests).

## Boundary conditions (by face-tag name, `bc_<tag> = wall | farfield | inflow | outflow`)
| kind | fluid | magnetic field |
|---|---|---|
| wall (planet surface, fluid/solid interface) | mirrored ghost (slip wall) | `E = 0` on the wall edges (flux through the wall frozen) |
| inflow (`farfield` with u.n_out < 0) | free-stream ghost state | `E = -u_inf x B_imf . dl` on its edges (a uniform stream is preserved exactly) |
| outflow (`farfield` with u.n_out >= 0) | zero-gradient ghost | free: one-sided face candidates, B evolves from the interior |
Solid cells (`n_solid` layers) are frozen: no flow, B frozen (E = 0).  The resistive solid shell is not implemented yet.
Defaults: tag `inner` = wall, `outer` = farfield; others `bc_default` (outflow).

## Output / restart
- ParaView: `DATA/sol_NNNNNNN_rRRRRR.vtu` + `.pvtu`, collected in `DATA/solution.pvd` (rho, p, u per species, B induced/total, J, u+).
- `DATA/fields_final.unsb`: global fields in the global numbering (`U_H`, `U_Na`, `B`, `rho`, `p`, `B_cell`, `J_cell`, `_time`).  It is the input of
  `mpcns_pre task=refine` (coarse solution) and of the solver itself (`init_fields=`; time is restored from `_time`).  The topology hash is stored as 0
  (the solver cannot hash the global complex); the refinement skips the check in that case.

## Verification (tests/test_solver.cpp, np = 1 and 4)
| test | result |
|---|---|
| A  uniform stream on a hanging-node polyhedral mesh, all faces inflow, 6 steps | change 1e-16 |
| B  Gaussian density wave, 3 resolutions | 2nd order: orders 2.11, 1.99 (L2 error 3.9e-4 at 160 cells); first order: 0.5–0.65 |
| C  Sod shock tube, 200 cells | L1 density error 3.5e-3, no overshoot |
| D  Alfven packet (B0 = 1, v_A = 1, t = 3) | peak 0.0249 at x = 2.975 (exact 0.0250 at 3.000); B_y = -u_y, div B ≤ 1e-13 |
| E  Mercury-like shell (dipole, wall, far field, solid core), 8 steps, np = 4 vs serial | max U difference 4e-16, Phi 3e-14, div B 4.5e-16 |
Pipeline (`tools/run_pipeline_demo.sh`): coarse run to t = 6 → refine from the coarse J and rho (1728 → 15322 cells, 1494 polyhedra, div B after the B transfer 2.5e-14)
→ run on the refined mesh to t = 8 (2452 steps, 8 ranks, div B 2.8e-14).

## Known limits of this version
- Physical fidelity is not validated against the legacy Mercury results (no Mercury benchmark was run); only the numerics above are verified.
- The cell-centred B is reconstructed with an area-weighted least-squares fit each stage (not cached); no performance tuning was done.
- Rusanov (first-order dissipation) with only the plasma velocity and fast speed in the EMF dissipation; no entropy fix.
- Hall limits (h² time step, sub-cycling) belong to the Hall work (docs/HallTodo.md).

Steady-state (pseudo-time) mode, local time stepping and the multigrid correction: see `docs/Multigrid.md` (`steady=1`, `mg_levels`, ...).
