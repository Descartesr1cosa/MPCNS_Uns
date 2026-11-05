# Mercury run: state of the work (saved mid-way — read this first when resuming)

Goal: the reference case `Run/CASE/setup` on a coarse mesh, to nondimensional time 25, on 28–56 ranks (limit: 84 processes),
without the Hall contribution to the EMF, with the fluid/solid load imbalance taken into account.

## Done
| item | where |
|---|---|
| Mesh: cubed-sphere solid shell r0..r1 (0.8..1.0, uniform layers) + fluid out to the eccentric ellipsoid x[-12,4], y,z[-10,10] (rays from the origin, geometric radial growth) | `CubedSphereSpec::{r_surf,ellipsoid,ell_lo,ell_hi}`, `Preprocess/INPUT/mercury_mesh.txt` (use `cs_r_stretch=1.07`: 1.15 gives 0.016-thick first cells on the sunward side) |
| Reference parameters in our format | `MercuryUns/INPUT/mercury.txt` (all values of `0_General_Setting … 3_MHD_Control`; Hall keys ignored) |
| B_add = IMF + dipole (SI moment/centre/softening converted like the legacy `DipoleField`), B_induced starts at 0; far-field inflow edges carry E = -u_inf x B_imf | `initial/`, `Config` |
| Na+ seed state (rho = 1e-6 rho_H, at rest, p = p_H 1e-6/23), exosphere (Exner et al. 2020), photo-ionisation rate, shadow, drag, thermal relaxation, thermal drift | `physics/Exosphere.*`, `fluid/FluidRHS.cpp::add_sources` |
| Ambipolar force on the ions (-w_s grad p_e) and ambipolar edge EMF (-pressure_coef grad p_e / n_e,eff, nodal pressure, d0 of a potential) | `physics/Ambipolar.*` |
| Mercury interior resistivity: implicit backward-Euler CT substep, `(star1/eta + dt D^T S2 D) E = D^T S2 Phi`, own parallel Jacobi-CG on the edge unknowns (no PETSc), window 0.84/1.04/0.02, eta_nd = eta_max/(mu0 U L) | `induction/Resistive.*` |
| Walls: H+ slip wall, Na+ absorbing (floor-state ghost); fluid/solid interface = wall, E = 0 on wall edges (the resistive substep then supplies E = eta J there); solid cells frozen, edges of solid-only cells skipped | `fluid/Reconstruction.cpp::ghost_state`, `boundary/` |
| Time step: legacy definition `CFL V / sum(|u.n| + c_fast) A`, magnetic wave speed from the **total** density (floor 1e-2), per-species sound speed | `fluid/FluidRHS.cpp::local_time_step` |
| Energy source with the mid-point velocity (`+ dt |S|^2 / 2 rho`) so the Lorentz force does not drain the thermal energy in low-beta cells | `add_sources` |
| Exact integration of the induce (ion gyration / pick-up) term by Strang splitting; conserves momentum and total kinetic energy | `physics/Gyration.*` |
| Timing bookkeeping for the load-balance calibration: `Solver::profile()` (wall, communication, owned fluid/solid cells, resistive edges, CG iterations), `Exchanger::comm_seconds()` | `driver/Solver.*` |
| `div B` is measured against the total flux (induced + background) | `diagnostics/DivB.*` |

All 14 unit/regression tests pass (`ctest`), including the earlier solver tests (the IMF now lives in the background field).

## Result with the first-order scheme (second_order=0)
Mesh `mercury_mesh.txt` (24192 cells, `cs_r_stretch=1.07`), 28 ranks, parameters of `mercury.txt`, `second_order=0`:
**reaches t = 25 stably** in 2 min 23 s (12190 steps, dt ~ 2e-3, divB <= 7.5e-14, min rho_H 0.033, min p_H 2.6e-3, Na+ mass 3.44 and growing,
H mass 3314 -> 3400). Command: add `second_order=0 output_step=5000` to the resume command below. Not yet analysed: the subsolar magnetopause
stand-off (the sampled x-axis profile still shows slow dense plasma, |u| < 0.3, down to the surface — the dayside pile-up / Na+ loading looks
strong and has not been checked against the legacy result), load balance, a comparison with MPCNS_Mercury output.

## Comparison with the legacy MPCNS_Mercury restart (`Run/../tmp/DATA`, step 764000, t = 24.51, 56 blocks, 524880 fluid cells)
Tools: `tools/legacy_compare/{rd.py,ref.py,ours.py}` (read the legacy `MPCNSRST` files + `.grd` grids; our `fields_final.unsb`). Ours: first order,
20736 fluid cells (25x fewer), t = 25.0; our output is kept in `tmp/uns_run1_first_order/`.
| quantity | legacy | ours (first order, coarse) |
|---|---|---|
| fluid volume / mean rho_H | 3640 / 1.035 | 3332 / 1.027 |
| mass H+ / Na+ | 3767 / 6.62 | 3421 / 3.46 |
| min rho_H / min p_H | 0.0061 / 0.0020 | 0.043 / 0.0049 |
| dayside axis: peak rho_H (x of peak) | 2.95 (1.6-1.8) | 2.38 (1.45-1.6) |
| dayside axis: u_x at x=1.8-2.2 / 1.6-1.8 | -0.81 / -0.25 | -0.79 / -0.42 |
| dayside axis: rho_Na at x = 1.0-1.1 / 1.2-1.3 | 0.77 / 0.43 | 0.10 / 0.09 |
| tail axis (x -2..-1.2): rho_H, p_H | 0.041, 0.054 | 0.083, 0.057 |
| tail axis (x -4..-3): rho_H, |u| | 0.060, 0.62 | 0.115, 0.84 |
The bow-shock pile-up, the pressure level and the deceleration towards the planet agree qualitatively and roughly in magnitude (the shock is
smeared over a thicker layer on our 25x coarser mesh); the tail density is ~2x higher and the wake less evacuated; the near-planet Na+ is 5-8x
lower and the total Na+ mass half (the exosphere scale heights 0.04-0.3 R_M are not resolved by our first fluid layers of 0.05 R_M, the hot-spot
sources are smeared, Na+ leaves through the absorbing wall).

## Load-balance profile (done)
`MercuryUns` writes `<output_dir>/load_profile.json` (full record) and `load_profile.txt` at the end of every run (`write_profile=1`):
per rank compute time/step (wall - communication/wait), wait time and fraction, owned fluid / solid cells, ghost ("virtual") cells and ghost/owned ratio,
owned edges, resistive edges, CG iterations; the least-squares fit `compute ~ a n_fluid + b n_solid` (solid/fluid ratio, r2), the wait per ghost cell,
and a `recommended_cost_model` (cost_fluid=1, cost_solid, comm_fluid, comm_solid). Use it for the refined mesh:
`mpcns_pre setup.txt mesh.txt cost_profile=<run>/load_profile.txt proc_num=N ...` (the .txt holds the `cost_*` keys; the partition report now also
lists the owned fluid / solid cells of each rank).
First measurements (24192 cells, 28 ranks, ~860 cells/rank, first order, 3 time units): compute max/mean = 1.11, **wait fraction 53-64 %**,
fitted solid cost 0.26-0.50 of a fluid cell (r2 0.2-0.4: rank-to-rank compute differences are small, so the fit is only a rough guide),
the repartition with the fitted weights did not change the balance measurably (1.11 -> 1.13, noise level). What matters more here:
- the implicit resistive CG: ~39 iterations/step with rtol 1e-6, **40 % of the run time** (19.8 s vs 11.6 s without it); `resistive_rtol=5e-3`
  (the reference value, now the default in mercury.txt) gives 15.7 s;
- at ~860 cells per rank the run is communication/latency bound (53 % wait even without the CG): fewer ranks (or a denser mesh) will be as fast per step.

## Refinement from the flow field and speed test (first order)
Coarse solution: `tmp/uns_run1_first_order/DATA/fields_final.unsb` (t = 25, 24192 cells). Refinement (`mpcns_pre task=refine`, max level 2, solid not refined, spherical snap):
manual shells `1.0-1.12 -> level 2`, `1.12-1.35 -> level 1`, plus indicators `h|J|/|B| > 1 -> 1` and `h|grad rho|/rho > 0.5 -> 1`
-> **124999 cells** (19359 / 29096 / 76544 at levels 0 / 1 / 2, 9743 polyhedra), first fluid layer 0.05 -> 0.0125 R_M, B transfer div B 3e-12.
The cell fields are prolonged piecewise-constant by default (`refine_prolong_linear=0`): the limited linear prolongation of the conservative
variables produced states that blew up in the first step (negative pressure / dt -> 1e-28). Command:
```
mpcns_pre setup.txt task=refine coarse_dir=<dir with global.unsb> coarse_fields=<fields_final.unsb> output_dir=<fine> proc_num=56 \
   refine_max_level=2 refine_solid_max_level=0 refine_spherical_snap=1 "refine_region_1=shell 1.0 1.12 2" "refine_region_2=shell 1.12 1.35 1" \
   "refine_indicator_1=current 1.0 1" "refine_indicator_2=density 0.5 1" [cost_profile=<coarse run>/load_profile.txt]
mpirun -np 56 MercuryUns setup.txt mercury.txt geometry_dir=<fine> init_fields=<fine>/fine_fields.unsb second_order=0 ...
```
Speed (600 steps, first order, dt ~ 3.8e-4, wall/step = slowest rank compute + wait):
| ranks | cells/rank | weights | wall/step | compute max/mean | wait |
|---|---|---|---|---|---|
| 28 | 4464 | default | 36.8 ms | 1.10 | 52 % |
| 28 | 4464 | coarse-run profile | 35.0 ms | 1.14 | 49 % |
| 56 | 2232 | default | 16.3 ms | 1.15 | 43 % |
| 56 | 2232 | coarse-run profile | 16.2 ms | 1.12 | 43 % |
| 84 | 1488 | default | 13.9 ms | 1.20 | 53 % |
=> 56 ranks is the sensible choice (84 gains only 15 %). Estimated cost of a run to t = 25 from the coarse state at this dt: 66000 steps x 16 ms ~ 18 min.
The coarse-run profile changes the step time by 0-5 % (within noise): at 2000-4000 cells per rank the balance is already ~1.1 and waiting
(latency, the resistive CG collectives) is what remains. Not run: the long refined run, second order on the refined mesh.

## NOT done yet
1. **The stability problem below** — no run has reached t = 25; the longest stable runs reach t ≈ 0.5–1.5 with features disabled.
2. Load balance: the profile and `cost_profile=` are in place (see above); still to do: judge it on the refined mesh (more cells per rank), and
   possibly restrict the CG collectives to the ranks that own resistive edges.
3. The production run (28–56 ranks, t = 25), its output, a magnetopause stand-off check (expected ≈ 1.4–1.5 R_M subsolar).
4. Unit tests for `Resistive`, `Gyration`, `Ambipolar`, `Exosphere`; update of `docs/Solver.md` (it still describes the first version).
5. Hall (EMF contribution) is intentionally not included.

## Open problem: instability of the full model at t ≈ 0.2–0.5 (mesh 24192 cells, 28 ranks)
Symptom: dt collapses to ~1e-10 within ~5 steps; the first sign is an anomalous Na+ pressure (e.g. p_Na ≈ 2e2 at rho_Na ≈ 2e-2,
T ≈ 1e4, one step before u_Na jumps) in first-fluid-layer cells (r ≈ 1.03–1.05), at several places near the cube-corner / panel-diagonal
directions and in the polar caps (B up to 600, rho_H ≈ 5, rho_Na up to 30 from the photo-ionisation hot spots).

Experiments (`max_Time` 0.4–2, CFL = legacy definition, 28 ranks; "ok" = reached the end):
| variation | result |
|---|---|
| baseline CFL 1.0, 0.5, 0.25, 0.12 (with exact gyration splitting) | fails at t ≈ 0.22 for all of them |
| `second_order=0` | ok (t=0.6) |
| `induce_scale=0` (no induce term) | ok to t=1.0 (min rho_H 0.02) |
| `chemistry=0` (no Na+ production, drag, relaxation) | ok to t=1.5 |
| `lorentz_force=0` | ok to t=1.5 (Na+ mass grows to 0.42) |
| `is_Ambipolar_Efield=0`, `is_Mercury_resistance=0`, `cuct=0`, gyro-frequency dt cap | still fail |
| old explicit induce term with CFL 0.25 | survives to t ≈ 0.45, then fails |
| `induce_scale=0.1` | ok to t=1.0 |
So the failure needs *all of*: Na+ production, the induce/Lorentz coupling and the second-order reconstruction. It is not a time-step
restriction (no CFL helps) and not a parallel bug (np = 4 and 28 fail too, differently). np-independence of the new physics has not been
re-tested over long runs.

Ideas, in order of cost:
- Na+ is a trace species with T and u set by sources: limit the reconstruction of the Na primitives (first order for rho_Na < small, or
  limit the temperature / use a pressure-positivity-preserving slope limiter), floor the Na temperature to ~T_n0.
- Check the energy source terms of Na at tiny density: `a5 n_s b1` and `a6 q T_n0` use the *true* n_s, the Lorentz mid-point term divides by rho_Na.
- Check the induce rotation for cells where rho_H or rho_Na is at the floor (chi from floored densities, f_H, f_Na).
- Test the pieces on a box (uniform B, Na source in a patch) to reproduce it without the planet; add a Na-only positivity test.
- Compare with the legacy scheme's handling: strong no-slip wall (first cell velocity zeroed), Na floor ghost rho = 1e-12, rho floor 2.5e-3 for H at the wall.

## How to reproduce / resume
```
export PATH=/workspace/software/openmpi-5.0.11/bin:$PATH
cmake -S . -B build && cmake --build build -j8 && (cd build && ctest)
W=/tmp/mercury_work; mkdir -p $W/run
build/Preprocess/mpcns_pre Preprocess/INPUT/setup.txt Preprocess/INPUT/mercury_mesh.txt output_dir=$W/geometry proc_num=28 cs_r_stretch=1.07 write_vtk=0
cd $W/run && mpirun -np 28 .../build/MercuryUns/MercuryUns MercuryUns/INPUT/setup.txt MercuryUns/INPUT/mercury.txt \
     geometry_dir=$W/geometry output_dir=$W/run/DATA max_Time=0.6 log_every=10        # fails near t = 0.22 (baseline)
# dumps around the failure:  output_step=1 output_start=55 max_Nstep=66 write_vtk=1  then
tools/inspect_vtu_extremes.py 0000060 $W/run/DATA
```
Run counts used so far: up to 28 MPI ranks, builds with -j8 (limit 84).
