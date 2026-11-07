# Mercury run: state of the work (read this first when resuming)

Last checked against the code: 2026-11-07 (second-order failure diagnosed, see "Diagnosis") (all 23 ctest tests pass). Sources moved after the split into `Common/` (model independent, namespace `mhd`) and
`MercuryUns/`: fluid fluxes/reconstruction/ghost states are in `Common/fluid/{FluxRHS,Reconstruction}.cpp`, the Mercury sources and time step in
`MercuryUns/fluid/MercurySources.cpp`, the Mercury-only induction pieces (`Resistive`) in `MercuryUns/induction/`, the generic implicit induction in
`Common/induction/ImplicitInduction.*`. Older rows below that name `fluid/FluidRHS.cpp` mean these files.
**Current default: `time_order=1` (forward Euler + semi-implicit induction).  Both `second_order=0` and `second_order=1` with `recon_rel_limit=0.25` (set in `mercury.txt`) reach t = 25 on the coarse mesh;
second order without the relative limiter fails at t ≈ 0.155 (see "Open problem" / "Diagnosis").**

Goal: the reference case `Run/CASE/setup` on a coarse mesh, to nondimensional time 25, on 28–56 ranks (limit: 84 processes),
without the Hall contribution to the EMF, with the fluid/solid load imbalance taken into account.

## Done
| item | where |
|---|---|
| Mesh: cubed-sphere solid shell r0..r1 (0.8..1.0, uniform layers) + fluid out to the eccentric ellipsoid x[-12,4], y,z[-10,10] (rays from the origin, geometric radial growth) | `CubedSphereSpec::{r_surf,ellipsoid,ell_lo,ell_hi}`, `MercuryUns/INPUT/mercury_mesh.txt` (use `cs_r_stretch=1.07`: 1.15 gives 0.016-thick first cells on the sunward side) |
| Reference parameters in our format | `MercuryUns/INPUT/mercury.txt` (all values of `0_General_Setting … 3_MHD_Control`; the Hall switch/taper are used by the semi-implicit induction below, Hall sub-cycling keys are ignored) |
| B_add = IMF + dipole (SI moment/centre/softening converted like the legacy `DipoleField`), B_induced starts at 0; far-field inflow edges carry E = -u_inf x B_imf | `initial/`, `Config` |
| Na+ seed state (rho = 1e-6 rho_H, at rest, p = p_H 1e-6/23), exosphere (Exner et al. 2020), photo-ionisation rate, shadow, drag, thermal relaxation, thermal drift | `MercuryUns/physics/Exosphere.*`, `MercuryUns/fluid/MercurySources.cpp::add_sources` |
| Ambipolar force on the ions (-w_s grad p_e) and ambipolar edge EMF (-pressure_coef grad p_e / n_e,eff, nodal pressure, d0 of a potential) | `physics/Ambipolar.*` |
| Mercury interior resistivity: implicit backward-Euler CT substep, `(star1/eta + dt D^T S2 D) E = D^T S2 Phi`, own parallel Jacobi-CG on the edge unknowns (no PETSc), window 0.84/1.04/0.02, eta_nd = eta_max/(mu0 U L) | `induction/Resistive.*` |
| Walls: H+ slip wall, Na+ absorbing (floor-state ghost); fluid/solid interface = wall, E = 0 on wall edges (the resistive substep then supplies E = eta J there); solid cells frozen, edges of solid-only cells skipped | `Common/fluid/Reconstruction.cpp::ghost_state`, `Common/boundary/` |
| Time step: legacy definition `CFL V / sum(|u.n| + c_fast) A`, magnetic wave speed from the **total** density (floor 1e-2), per-species sound speed | `MercuryUns/fluid/MercurySources.cpp::local_time_step` |
| Energy source with the mid-point velocity (`+ dt |S|^2 / 2 rho`) so the Lorentz force does not drain the thermal energy in low-beta cells | `add_sources` |
| Exact integration of the induce (ion gyration / pick-up) term by Strang splitting; conserves momentum and total kinetic energy | `physics/Gyration.*` |
| Timing bookkeeping for the load-balance calibration: `Solver::profile()` (wall, communication, owned fluid/solid cells, resistive edges, CG iterations), `Exchanger::comm_seconds()` | `driver/Solver.*` |
| `div B` is measured against the total flux (induced + background) | `diagnostics/DivB.*` |

All 23 ctest tests of the repository pass (21 library/Mercury, 2 Mars), including the earlier solver tests (the IMF now lives in the background field).
There are still no dedicated unit tests for `Resistive`, `Gyration`, `Ambipolar`, `Exosphere` (they are only exercised by `test_solver`/`test_implicit_induction` and the full runs).

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

## Refined run t = 25 -> 30 and speed against the legacy program
Run: refined mesh above (124999 cells), first order, 56 ranks, restart time carried through the refinement (`_time`), `cost_profile` of the coarse run.
Result in `tmp/uns_run2_refined_first_order/`. Wall time **7 min 16 s** for 30127 steps (dt 3.8e-4 -> 1.4e-4 as the cells near the planet fill with Na+), i.e.
**14.5 ms/step, 87 s per unit of nondimensional time**, divB <= 4.2e-13, stable.
Legacy logs (`tmp/DATA/M2*.o*`, 56 ranks, 524880 cells, dt = 3.2e-5 constant):
| | legacy MPCNS_Mercury | MercuryUns (refined, first order) |
|---|---|---|
| cells (fluid) | 524880 (0.525 MDOF) | 124999 (115495 fluid) |
| dt | 3.2e-5 | 3.8e-4 .. 1.4e-4 |
| wall per step | 0.20 s (run 1, t 0-60), 0.114 s (continuations) | 14.5 ms |
| wall per cell-step | 3.9e-7 s (run 1) / 2.2e-7 s (continuations) | 1.16e-7 s |
| wall per unit time | 6900 s (run 1) / 3550 s (continuations) | 87 s |
Per step the new code is 8-14x faster, per cell-step 2-3x faster, and per unit of simulated time 40-80x faster (4x fewer cells and a 4-10x larger dt).
Caveats: different machines (legacy: InfiniBand cluster, Intel; here: one 256-core node, shared memory MPI), the legacy logs include its I/O, implicit
resistive PETSc solves and Hall subcycling settings, our scheme is first order with a coarser mesh (the legacy flux is also first-order Rusanov). The
per-cell figure is the fairest one; the per-time figure mostly reflects the mesh and dt.
State at t = 30 against the legacy t = 24.5 (near the axis, cylinder radius 0.3): dayside x 1.0-1.1 rho_H 0.36 / 0.28, rho_Na 0.23 / 0.87;
x 1.3-1.6 rho_H 2.12 / 2.32; x 1.6-1.8 u 0.37 / 0.26; tail x -2..-1.2 rho_H 0.052 / 0.041, x -4..-3 rho_H 0.085 / 0.057 (ours / legacy);
Na+ mass 3.94 / 6.62, min rho_H 0.0102 / 0.0061. The refinement moved everything toward the legacy values (coarse run: Na near the planet 0.10, tail
rho_H 0.083, min rho_H 0.043) but near-planet Na+ is still 4x lower.

## Semi-implicit induction: ambipolar, Hall and Mercury resistivity (`time_order = 1`)
Method (`induction/ImplicitInduction.*`, test `tests/test_implicit_induction.cpp`):
- unknown = the total edge EMF `E` (1-cochain); `Phi^{n+1} = Phi^n - dt d1 E` is the only magnetic update, so **div B is preserved for any E, i.e. for any
  linear-solver tolerance** (measured 1e-14 .. 4e-13 with rtol 1e-3);
- `E = E_expl + F(Phi^{n+1})`, `E_expl` = ideal (UCT, arbitrary valence) + ambipolar EMF at t^n, boundary values (wall E = 0 before the implicit part, inflow E = -u x B_imf imposed
  as Dirichlet, not an unknown); `F = M d1^T S2 Phi + H(Phi)`;
- resistive part `M = eta/star1` (Mercury window 0.84/1.04/0.02, eta_nd = eta_max/(mu0 U L)), implicit in `Phi^{n+1}`;
- Hall part `H(Phi)_e = sum_c w_{e,c} taper alpha_c (J_c[Phi] x B0_c).t_e`, `alpha_c = hall_coef/n_e,eff`, radial taper r_min..r_max = 1.0..1.012 (reference), `B0` = total cell field at
  t^n frozen (linearised Hall), `J_c[Phi]` the Gauss curl of the cell field reconstructed from Phi; edge sum with the dual-face weights through the owner reduction (arbitrary valence);
- linear system `(I + dt F d1) E = E_expl + F(Phi^n)`: nonsymmetric (Hall ~ skew, resistive SPD, spectrum in the right half plane) -> restarted GMRES(20), Jacobi right
  preconditioner, matrix-free (one edge exchange per matvec); np-independent;
- fluid: forward Euler; ion gyration by exact Strang-split rotation; Na+ photo-ionisation / drag / relaxation stay explicit: their rates (drag a4 vst = 3e-4, relaxation
  a5.. ~1e-4, production 7 per unit time into a trace seed density) are 1e-4 .. 1e-1 of 1/dt (dt ~ 1e-4), i.e. not stiff — no implicit treatment added.
Checks: Gaussian whistler packet (`omega = alpha B0 k^2`, `E_H = alpha J x B0`) against the exact Schroedinger-type solution: L2 error 0.19 (h=0.08) -> 0.072 (h=0.04), time error
decreasing with dt; stable with dt = 77x the explicit whistler limit (energy only decreases); np = 1 and 4 identical to 1e-17.
Mercury case: sparse mesh 28 ranks, t = 3: Hall off/on both stable, GMRES 19 it/step at rtol 1e-3 (47 at 1e-6, solution identical to 5 digits), 20 s per 3 time units.
Refined mesh (125k cells, 56 ranks, t = 25 -> 26, `mercury.txt`): Hall off 3579 steps / 121 s, Hall on 3578 steps / 132 s (+9 %), 27 GMRES it/step,
33.6 / 36.6 ms per step. The explicit whistler limit printed by the solver (`dt_hall_explicit`) is 5e-5 against fluid steps of 1e-4 .. 2e-4, so an explicit Hall term would cost
2-4x more steps; the Hall term changes the state at t = 26 by only ~1 % (hall_coef = 1.3e-4: Hall is a small correction in these units).
Not done: second-order time integration of the induction, Hall sub-cycling (not needed), tuning of the preconditioner (the Jacobi preconditioner ignores the L coupling).

## Preconditioning of the implicit induction, dense mesh and stability test
Measured on the refined mesh of 125k cells (56 ranks, 1500 steps from t = 25, rtol 1e-3 unless noted):
| variant | GMRES it/step | wall/step | max rel. difference to the tight solution (rtol 1e-9) |
|---|---|---|---|
| Jacobi, cold start, rtol 1e-9 (reference) | 72.8 | 97.7 ms | - |
| Jacobi + **warm start** (previous step's implicit correction as initial guess) | 2.7 | 12.7 ms | U_H 4e-5, B flux 4e-7 |
| shifted block-Jacobi ILU(0) + warm start | 1.6 | 14.1 ms | U_H 5e-4, U_Na 2e-3, B flux 1e-5 |
On the sparse mesh (28 ranks) ILU(0) was 1.7x better than Jacobi in iterations (2.1 vs 3.6 it/step), on the refined mesh it is not (more setup, weaker for the
curl-curl part; a plain ILU(0) of `I + dt M L` is unstable for large dt, hence the diagonal shift 0.3). Conclusion: the warm start is what matters (19 -> 3.6 it/step on
the sparse mesh, 73 -> 2.7 on the refined one; the old explicit step cost 16 ms, the semi-implicit one with Hall now 12.7 ms); defaults: Jacobi, rtol 1e-4 (`induction_ilu=1` selects ILU).
Hall explicit limit for reference: `dt_hall_explicit` in the log (1e-4 .. 1e-5 here).

Dense mesh (`scripts/03_refine.sh`): 574665 cells (14114 / 46647 / 237232 / 276672 at levels 0-3, 45924 polyhedra), first fluid layer 0.0063 R_M (8x finer than the coarse
mesh), B transfer div B 1.1e-11, partition imbalance 1.02 (cost profile of the coarse run). Short run from the coarse state at t = 25, 56 ranks, 3000 steps:
**stable** (4 min 22 s): divB 1.6e-12 constant, min rho_H 0.043 -> 0.028 (the cavity keeps deepening while Na+ accumulates; the legacy minimum is 0.006), dt 1.5e-4 -> 8.5e-5
(still decreasing slowly as the near-planet cells fill up), 9 GMRES it/step, 87 ms/step, compute imbalance 1.12, wait 28 %. Estimated cost of 25 more time units:
~3e5 steps x 87 ms = 7 h. Scripts for the whole chain: `scripts/` (README there).

## NOT done yet (status 2026-11-07)
1. The second-order stability problem is fixed on the coarse mesh by `recon_rel_limit` (below); still to do: second order on the refined meshes, comparison with the legacy field.  Older statements "no run has reached t = 25" are obsolete.
2. Load balance: judged on the refined mesh (section above): imbalance 1.02–1.2, waiting (latency) dominates; the restriction of the resistive-CG collectives
   to the ranks that own resistive edges is only relevant for `time_order=2` and was not done.
3. A second-order production run (t = 25) and its comparison with the legacy run; a magnetopause stand-off check (expected ≈ 1.4–1.5 R_M subsolar).
4. Unit tests for `Resistive`, `Gyration`, `Ambipolar`, `Exosphere`.  (`docs/Solver.md` was brought up to date on 2026-11-07.)
5. Hall is **included** since the semi-implicit induction (section above); not done: Hall sub-cycling (not needed) and second-order time accuracy of the Hall/resistive terms.

## Open problem (SOLVED 2026-11-07 by `recon_rel_limit`, see "Diagnosis"): instability of the full second-order model at t ≈ 0.2–0.5 (mesh 24192 cells, 28 ranks)
Symptom: dt collapses to ~1e-10 within ~5 steps; the first sign is an anomalous Na+ pressure (e.g. p_Na ≈ 2e2 at rho_Na ≈ 2e-2,
T ≈ 1e4, one step before u_Na jumps) in first-fluid-layer cells (r ≈ 1.03–1.05), at several places near the cube-corner / panel-diagonal
directions and in the polar caps (B up to 600, rho_H ≈ 5, rho_Na up to 30 from the photo-ionisation hot spots).

Re-checked on 2026-11-07 with the current code (12 ranks, same 24192-cell mesh, `time_order=1` default of `mercury.txt`, `second_order=1`): dt collapses to 1e-25 at
t = 0.1554 (step ~55), `min p_H` -> 1e-12, `mass_Na` 0.21 -> 4e2; `second_order=0` runs to t = 0.6 with `mass_Na` 0.18.  The experiments in the table below
were done earlier, mostly with the SSP-RK2 path (`time_order=2`, explicit EMF); the failure time is the same order in both time integrators, so the cause is not the
time integration of the induction.
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

### Diagnosis (2026-11-07, 12 ranks, each run takes ~3 s; temporary debug hooks, reverted)
| experiment (`t_end` 0.6–3) | result |
|---|---|
| second order for H+ only (Na+ first order) | **stable** to t = 0.6 |
| second order for Na+ only (H+ first order) | fails at t = 0.126 → the instability lives in the Na+ reconstruction |
| count of floor events (`enforce_floors`) | Na+ **density becomes negative** already at step 5 (t = 0.019, e.g. rho = -1e-5 .. -1.7e-4 in first-layer cells, r ≈ 1.03) and the cell is reset to the floor state; 4762 events by t = 0.16; none with `second_order=0` |
| Na+ mass | identical to first order until t ≈ 0.08, then 3x larger at t = 0.154 (0.21 vs 0.066): the floor reset replaces a negative mass by +1e-10 and the momentum/energy by the floor state, i.e. it is a spurious source/sink, not physics |
| first anomaly in T_Na = p/rho | T_Na = 0.59 → 11 → 1.8e3 → 1e5 in 12 steps (t = 0.124 → 0.155) in dayside first-layer cells where H+ is compressed (rho_H ≈ 2) and B ≈ 150; correct value ≈ T_n0 = 6e-3 |
| extra *relative* limiter on the reconstruction of rho and p, face value within `[(1-θ) q_c, (1+θ) q_c]` (Na+ only) | θ = 0.5 / 0.25: Na+ stable (no Na+ floor events), run to t = 1 with mass_Na 0.29; θ = 1: still fails |
| the same for H+ and Na+ (θ = 0.5, 0.25) | stable to t = 3 (mass_Na 0.69 / 0.67, min rho_H 4e-3 / 8e-3); θ = 0.5 still has ~400 H+ *pressure* floor events (E ≈ ke to 4 digits, r ≈ 1.09), θ = 0.25 fewer |

**Cause.** The Barth–Jespersen limiter in `Stencil::limited_gradients` bounds the reconstructed face value by the *neighbour* extremes, not by the cell value.
Na+ is a trace species whose density jumps by 1e2–1e4 between neighbouring first-layer cells (seed 1e-6 against the photo-ionisation hot spots),
so a trace cell can be extrapolated to a face value tens to thousands of times its own mean.  The forward-Euler update of that cell then removes
`dt/V sum_f (rho_f u + s rho_f) A` with `rho_f >> rho_c`: the density (and the thermal energy, since `p = (gamma-1)(E - m²/2rho)` is a difference of two nearly equal numbers)
goes negative, the floor reset injects mass and kills momentum, and the cell next to it sees a heated, over-dense neighbour (T_Na up by 1e4–1e8) → `dt` → 1e-25.
Positivity of a MUSCL step needs `CFL * rho_f/rho_c <~ 1` (Zhang–Shu type condition); the CFL here is the *cell-wise* sum of wave speeds with `rho_f ≈ rho_c` assumed,
so lowering CFL (0.12 failed too) cannot repair a ratio of 1e2–1e4; the first-order scheme has `rho_f = rho_c` and is positivity preserving under CFL ≤ 1.
It needs all of Na+ production (creates the large contrasts), the induce/Lorentz coupling (moves Na+ and H+ states apart in the strongly magnetised layer) and the second-order
reconstruction — exactly the table above.  Not a parallel bug (np = 4, 12, 28 all fail) and not caused by the time integrator of the induction.
Secondary effect found on the way: low-beta H+ cells (B ≈ 150, p_H ≈ 1e-4 of the kinetic energy) lose their pressure to round-off in `E - ke`; first-order runs do not show it.

**Fix implemented (2026-11-07): `recon_rel_limit` = θ** (`Stencil::set_rel_limit`, `Common/fluid/Reconstruction.cpp`; default 0 = off, so Mars is unchanged; `mercury.txt` sets 0.25).
After Barth–Jespersen, `psi <- min(psi, θ |q_c| / |g.dx_f|)` for rho and p over all faces.  Result, coarse mesh (24192 cells), 12 ranks, `second_order=1`, `time_order=1`, θ = 0.25:
**reaches t = 25 stably** (28247 steps, 3 min 56 s, divB 6e-14, min rho_H 4.7e-3, min p_H 1.3e-3, mass_H 3404, mass_Na 3.87; dt 6.9e-4 at t = 24.8).
For comparison, first order: min rho_H 0.033, mass_Na 3.44; legacy (t = 24.5): min rho_H 0.0061, mass_Na 6.62 — the second-order cavity/minimum density is much closer to the legacy value.
θ = 0.5 also ran to t = 3 (more H+ pressure floors); θ = 1 and off fail.  All 23 ctest pass.  Not yet done: field comparison with the legacy run for the second-order state,
the H+ internal-energy/low-beta fix (item 3 below), a refined-mesh second-order run, a unit test of the limiter.
**Original recommendation (item 1 is what was implemented):**
1. Positivity-preserving limiter on the face values of rho and p: `psi_c <- min(psi_c, theta |q_c| / |g.dx_f|)` over the faces, theta ≈ 0.25–0.5 (the experiment above), or the exact Zhang–Shu bound from the CFL;
2. keep p and T consistent (reconstruct T or limit rho, p with one common psi) and fall back to first order in a cell if the face states violate `rho_f > 0.2 rho_c` or `p_f > 0.2 p_c` after limiting;
3. for the H+ low-beta cells: evolve the internal energy (or use the entropy/pressure fix of the Mars branch) instead of recovering `p` from `E - ke`.
Then re-run the second-order case to t = 25 and compare with the first-order result and the legacy run.

Ideas, in order of cost (earlier list; items 1 and 2 were superseded by the diagnosis above):
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
build/Preprocess/mpcns_pre Preprocess/INPUT/setup.txt MercuryUns/INPUT/mercury_mesh.txt output_dir=$W/geometry proc_num=12 cs_r_stretch=1.07 write_vtk=0
R=$PWD; cd $W/run && mpirun -np 12 $R/build/MercuryUns/MercuryUns $R/MercuryUns/INPUT/setup.txt $R/MercuryUns/INPUT/mercury.txt \
     geometry_dir=$W/geometry output_dir=$W/run/DATA max_Time=0.6 log_every=10        # second order: fails near t = 0.155 (3 s wall); add second_order=0 for the stable run
# dumps around the failure:  output_step=1 output_start=55 max_Nstep=66 write_vtk=1  then
tools/inspect_vtu_extremes.py 0000060 $W/run/DATA
```
Run counts used so far: up to 28 MPI ranks, builds with -j8 (limit 84).
