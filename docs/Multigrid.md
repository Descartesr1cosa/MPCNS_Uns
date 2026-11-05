# Multigrid / pseudo-time acceleration of MercuryUns (steady-state mode)

Status (2026-11-05): implemented, tested (20 ctest tests), measured on the coarse (24192 cells) and a refined (125 k cells) Mercury mesh; every experiment is recorded, including the failed ones (section 9).

**Summary of what was found**
* A steady-state mode (`steady=1`) with *local pseudo-time stepping* plus an agglomeration *FAS multigrid* correction of the fluid was added; the time-accurate solver is unchanged.
* With the induction step capped at R = 32 × the smallest cell step the iteration converges in ≈ 1150 iterations (≈ 21× fewer than the 24000 steps of the time-accurate run; the multigrid halves the iterations of local stepping alone and lowers `res_H` by two orders at equal iteration count),
  **but** the converged state is 15–70 % away from the true steady state: with a per-edge step the fixed point of constrained transport is $d_1(TE)=0$, not $d_1E=0$ (section 2, table 6.5).
* R = 4 gives a few-percent state with ≈ 8× fewer iterations / 6× less CPU; there the convergence is limited by **B** and the multigrid (which acts on the fluid only) brings nothing (6.6).
* On the refined mesh the multigrid removes the stagnation of the local-step iteration (res_H 3e-4 instead of 1.5e-3 after 3000 iterations, 6.3).
* Open problem / next step: accelerate B (B hierarchy from the nested octree, or implicit ideal induction), section 8.

## 1. What is being accelerated, and why this design

The Mercury production run advances in *physical* time with one global step `dt = min_c dt_c`.  On the refined meshes the smallest cells
(first fluid layer 0.006–0.0125 R_M) in the strongly magnetised polar caps (|B| up to 600 in units of the IMF) give `dt ~ 1e-4`, whereas cells in the
far field could take `dt ~ 0.3`: a spread of 10^2.5.  To reach the quasi-steady magnetosphere (t = 25 and beyond; the legacy code needed 7·10^5 steps) the
global step wastes almost all of its work on cells that could move 300× faster.  (On the coarse mesh the global-step run does reach the steady state,
but only at t ≈ 45–75, i.e. 22000–37000 steps; see 6.5.)  The *steady-state mode* (`steady=1`) therefore iterates a pseudo-time,
and a multigrid correction removes the smooth part of the error that the (explicit, local) smoother damps slowly.

Three ingredients, each switchable and documented below:

| # | ingredient | switch | what it does |
|---|---|---|---|
| 1 | local pseudo-time stepping (LTS) | `steady=1` | every fluid cell uses its own CFL step; the induction uses one step per edge |
| 2 | step cap for the induction | `lts_b_ratio=R` | edge step ≤ R · (smallest cell step); needed for the stability of CT with a variable step (section 3) |
| 3 | FAS multigrid for the fluid | `mg_levels=L` | agglomeration hierarchy of clusters, first-order coarse operator, B frozen on the coarse levels |

Time-accurate runs (`steady=0`) are untouched (regression-checked bit-wise: the global-step path of `step_first_order` is the old code plus residual bookkeeping).

## 2. Local pseudo-time stepping

The semi-discrete system (two ion species, constrained transport) is
$$\partial_t U_s = R_s(U,\Phi) ,\qquad \partial_t \Phi = -d_1 E(U,\Phi) ,\qquad d_2\Phi = 0 .$$
Steady states are the zeros of the right-hand side.  Replacing $dt$ by a cell-dependent $\tau_c$ does not change the zeros of the fluid part:
$$U_c^{m+1} = U_c^m + \tau_c\,R_c(U^m),\qquad \tau_c = \mathrm{CFL}\;\frac{V_c}{\sum_{f\in\partial c}(|u_s\!\cdot n_f| + c_{\rm fast})A_f}\ \text{(species minimum)},$$
and the iteration is still a convex combination (stable) cell by cell because each cell satisfies its own CFL bound.

**Constrained transport.**  The magnetic update has to keep $d_2\Phi=0$ exactly.  With a face-dependent step $\Phi_f -= \tau_f (d_1E)_f$ the divergence
is lost.  Using a per-edge step instead keeps it: the update is
$$\Phi^{m+1} = \Phi^m - d_1(T E),\qquad T=\mathrm{diag}(\tau_e),\qquad \tau_e=\min_{c\ni e}\tau_c ,$$
and $d_2 d_1 = 0$ holds for any 1-cochain $TE$ (implemented in `faraday_update_edge_dt`).  The implicit part (resistive shell, Hall) enters through
$$(I + F\,d_1 T)\,E = E_{\rm expl} + F(\Phi^m),\qquad F(\Phi)=M d_1^{\!\top}\!S_2\Phi + H(\Phi),$$
solved by the warm-started GMRES with the Jacobi diagonal $1+\tau_e m_e d_{2,e}$ (`ImplicitInduction::solve(..., dt_edge)`).

**Where the fixed point is no longer exact.**  The fixed point of the B-update is $d_1(TE)=0$, not $d_1E=0$.  A steady ideal-MHD flow has
$E=-\nabla\varphi\neq0$ (curl free), so $d_1(TE)=0$ requires $\tau_e$ to be uniform wherever $E\neq0$.  With a strongly varying $\tau_e$ the LTS-CT
iteration therefore converges to a *different* state (error $\sim |\nabla\tau|\,|E|/\tau$) or does not converge at all.  Measured (section 6): with an
unlimited edge step the run develops a global instability (mass of H+ grows by a factor 15); the *true* residual
$\|\Phi^{m+1}-\Phi^m\|/(\tau_f\|\Phi\|)$ (not $d_1E$!) must be used for monitoring.  Capping the edge step at $R\,\tau_{\min}$ makes $\tau_e$ uniform ($=R\tau_{\min}$) in all
cells where the CFL step exceeds the cap, i.e. in most of the domain, so the error is confined to the transition region around the planet.
This is the compromise used here (`lts_b_ratio`, default 16 in the examples); the fluid keeps its full local step.

## 3. Residual definitions used in the logs
* `res_H`, `res_Na`: $\big(\sum_c V_c (\Delta\rho_c/\tau_c)^2\big)^{1/2}/\big(\sum_c V_c\rho_c^2\big)^{1/2}$ (rate of change of the mass density, from the explicit right-hand side).
* `res_B`: $\big(\sum_f(\Delta\Phi_f/\bar\tau_f)^2\big)^{1/2}/\big(\sum_f\Phi_{{\rm tot},f}^2\big)^{1/2}$ over the faces touching fluid cells, `res_Bsolid` the same for solid-only faces.
Both are rates per unit pseudo time, so global-step and local-step runs are comparable.

## 4. FAS multigrid for the fluid (`multigrid/Multigrid.*`)

### 4.1 Hierarchy (agglomeration, per rank, no new mesh)
Level 0 is the set of owned fluid cells.  Level $l+1$ is built from level $l$ by **heavy-face matching** on the cluster graph:
clusters $I,J$ of one rank are paired in order of the weight $w_{IJ}=\|\sum A_{IJ}\|$ (shared area vector); unmatched clusters join their most strongly connected neighbour;
`mg_agglom_passes` (default 2) pair-matching passes per level give clusters of ~4 members (20736 → 5151 → 1266 → 309 cells on the coarse Mercury mesh).
Heavy faces first means thin cells (first fluid layers, anisotropic) are merged across their large faces, i.e. along the thin direction (semi-coarsening).
Clusters never cross a rank boundary, so the hierarchy needs no repartitioning and is a pure function of the partition.

For a cluster $I$
$$V_I=\sum_{c\in I}V_c ,\qquad A_{IJ}=\sum_{f\in\partial I\cap\partial J}A_f\ (\text{outward from }I),\qquad A_{I,k}^{\rm bnd}=\sum_{f\in\partial I,\ \text{kind }k}A_f ,$$
and faces are merged only if their normals agree ($\cos>0.9$).  (The first version merged all boundary faces of the same kind: a corner cluster then
got one *diagonal* face of area $\sqrt3$ instead of $3$, the level-0 flux differed from the fine first-order right-hand side by 35 % and the 4-level cycle diverged;
found by `test_multigrid` B.)  Closure $\sum_J A_{IJ}+\sum_k A^{\rm bnd}_{I,k}=0$ holds to round-off.
Clusters owned by other ranks that touch an owned cluster ("remote clusters") are exchanged through the existing fine-cell `Exchanger`: every owned fine cell
carries the value of its cluster, the ghost fine cells hand the remote values back (no new MPI plan).

### 4.2 Coarse operator
Rusanov flux of the cluster averages (first order), the same source terms as on the fine level with the frozen restricted fields $\bar B,\bar J,\overline{\nabla p_e}$,
$\bar q_{\rm photo},\bar\nu$ and the same boundary ghost states (`ghost_state`: slip wall for H+, absorbing wall for Na+, free stream at inflow):
$$R_I(U)=-\frac1{V_I}\Big[\sum_J\Phi(U_I,U_J;n_{IJ})|A_{IJ}|+\sum_k\Phi(U_I,U^{\rm ghost}_k;n_k)|A_k^{\rm bnd}|\Big]+S_I(U_I;\bar B,\bar J,\dots).$$
One **smoothing step** on any level is the pseudo-time step of the fine solver (section 2) with a cluster-wise CFL step,
$$\mathcal S_l(U)=G_{\tau/2}\circ(\mathrm{id}+\tau R_l)\circ G_{\tau/2}(U),$$
$G_{\tau}$ = exact gyration of the relative ion velocity about $\bar B$ (`gyrate_n`).  Its pseudo rate is $\rho_l(U)=(\mathcal S_l(U)-U)/\tau$.

### 4.3 FAS cycle
With $\bar I$ = volume-weighted restriction (conservative: $\sum_I V_I\bar U_I=\sum_c V_cU_c$) and injection as prolongation:

1. fine state $u_0$ after the pseudo step, fine rate $\rho_0(u_0)$ evaluated on a copy (fluid part only: gyration, flux, sources; no induction);
2. level $l=1\ldots L-1$ (down): $u_l^{\rm in}=\bar I u_{l-1}$, restricted rate $r_l=\bar I\,\rho^{\rm tot}_{l-1}$, forcing
   $$g_l=r_l-\rho_l(u_l^{\rm in}),$$
   $\nu_1$ smoothing steps $u_l\leftarrow\mathcal S_l(u_l)+\tau\,g_l$ (the forcing is added *outside* the Strang splitting, so the fixed point is exact), then
   $\rho^{\rm tot}_l=\rho_l(u_l)+g_l$ and restriction to $l+1$; the coarsest level does `mg_coarse` steps;
3. up: $e_l=u_l-u_l^{\rm in}$, $u_{l-1}\mathrel{+}=\alpha\,e_l$ on the members, then $\nu_2$ post-smoothing steps (not on the fine level: the next fine step is the post-smoother);
4. the correction $e_1$ goes to the fine cells with `positivity_safeguard`: $\alpha\in\{1,\tfrac12,\tfrac14,\dots\}$ the largest value keeping $\rho\ge0.1\rho_{\rm old}$,
   $e_{\rm int}\ge0.1\,e_{\rm int,old}$, $\rho\le$ `mg_growth` $\rho_{\rm old}$ (and the same for the internal energy); the safeguard is cell-wise and is the reason the cycle never produced a negative pressure.

**Fixed point.** If $\rho_0(u_0^*)=0$ then $u_l=\bar I u_{l-1}^*$ makes $\rho_l+g_l=\bar I\rho_{l-1}^{\rm tot}=0$ on every level: the correction vanishes and the converged fine solution
is a fixed point of the cycle, independent of the quality of the coarse operator (`test_multigrid` C, D).

**B is frozen on the coarse levels.**  The induction (CT, implicit Hall/resistive solve) stays on the fine level; the coarse levels see $\bar B$, $\bar J$ as data.  The cycle therefore
accelerates the fluid; B converges at the rate set by its own pseudo step (section 2) and by the slow fluid modes it follows.  Measured consequence in section 6: $res_B$ is the last quantity to converge.

**Cost.** One cycle costs ~(2 + 1) smoothing sweeps on a hierarchy that shrinks by ~4 per level plus the exchanges (each coarse smoothing step performs one fine-cell exchange);
a cycle every `mg_cycle_every`=4 fine steps adds ≈ 0–10 % per step.

### 4.4 Parameters
| key | default | meaning |
|---|---|---|
| `steady` | 0 | 1: pseudo-time iteration (local steps) instead of the time-accurate global step |
| `lts_b_ratio` | 0 (none) | cap of the edge (induction) step in units of the smallest cell step; **use 16–32** |
| `lts_max_ratio` | 0 | cap of the cell step too (diagnostic) |
| `lts_solid_dt` | 1 | step of edges without fluid cells (resistive shell), in units of the largest cell step |
| `lts_smooth` | 0 | gradient limiter of the step field (tested, not helpful, see log) |
| `mg_levels` | 0 | number of levels incl. the fine one (0/1: off).  4 recommended, 5 diverged |
| `mg_cycle_every` | 1 | fine steps between corrections (4 recommended) |
| `mg_start` | 0 | first step with a correction (200 recommended: the initial transient is source driven) |
| `mg_pre`, `mg_post`, `mg_coarse` | 1, 1, 2 | smoothing steps; more smoothing *diverged* (extrapolation too long) |
| `mg_relax`, `mg_growth` | 1, 2 | damping and growth cap of the correction |
| `mg_agglom_passes`, `mg_min_clusters` | 2, 8 | agglomeration |
| `steady_tol` | 0 | stop when `res_H` and `res_B` are below it |
| `debug_res` | 0 | print the face with the largest B rate and non-finite coarse states |

## 5. Tests (all in `ctest`, np = 1 and 4; 20 tests pass)

| test | what is checked | result |
|---|---|---|
| `test_multigrid` A | hierarchy on a hanging-node mesh (2633 cells → 655 → 162 → 40): cluster volumes add up, every cluster surface closes, total area vector between two clusters is antisymmetric, restriction conserves $\sum V U$ | 0 / 0 / 0 (closure, antisymmetry, volume), conservation 5e-15 |
| `test_multigrid` B | level 0 (clusters = cells) flux right-hand side equals `fluid_flux_rhs` with zero gradients, np = 1, 4 (includes remote clusters) | 2.2e-15 (max rhs 5) — this test found the diagonal-boundary-face bug (35 % error) |
| `test_multigrid` C | uniform stream, zero rate → the cycle returns exactly zero (fixed point) | 0 |
| `test_multigrid` D | zero smoothing steps → zero correction for an arbitrary smooth state and rate (the forcing $g$ is consistent with the restricted rate) | 0 |
| `test_multigrid` E | positivity safeguard (decrease, growth cap) | alpha = 0.25 / 1 / 0.5 |
| `test_pseudo_time` A | `faraday_update_edge_dt`: $d_2B$ after 20 updates with random edge steps (factor 15 between edges) on a hanging-node mesh | 5e-16 |
| `test_pseudo_time` B | steady mode with `lts_max_ratio=1` reproduces the global-step first-order scheme (8 steps, two ion species, Lorentz force, CT) | 2.2e-16 |
| `test_pseudo_time` C | local steps (range ×4.6), 40 iterations: $d_2 B$ 1e-16, density positive | pass |
| old tests | `test_solver`, `test_ideal_emf`, `test_implicit_induction`, ... after the refactors (`add_sources_n`, `gyrate_n`, `cell_time_steps`, per-edge `ImplicitInduction::solve`) | unchanged, all pass |
| np independence of the converged state | coarse mesh, 2000 iterations, LTS + 4-level MG, 3 ranks vs 7 ranks (different partitions → different hierarchies) | mass_H 3408.760 vs 3408.768, mass_Na 3.12766 vs 3.12761, min rho_H identical (3.277e-2) |

The time-accurate path is unchanged: the coarse run to t = 25 on 12 ranks reproduces the earlier numbers (mass_H 3399.8, min rho_H ~0.04).

## 6. Results

All timings: this 256-core node while another 56-rank production run was using it (load ≈ 60–70), so the absolute wall times are noisy (±15 %), the iteration counts are exact.
Test processes: never more than 15 at a time (3-rank runs ×4, 7-rank runs ×2, ctest ≤ 7).
Metric: first logged iteration with `res_H < 1e-5` **and** `res_B < 1e-4` (`tools/conv_report.py`).  Logs: `tmp/mg_work/runs/<name>/run.log`.

### 6.1 Coarse mesh (24192 cells, 3 ranks, `second_order=0`, start from the analytic initial state, `lts_b_ratio=32`, log every 50)

| run | configuration | iterations to tolerance | wall to tolerance | final res_H / res_B (2000 it) |
|---|---|---|---|---|
| qA | LTS only | not reached (≈ 2300 by extrapolation) | ≈ 90 s | 5.1e-5 / 9.2e-5 |
| qB | LTS + MG, 3 levels, cycle every step | 1050 | 53 s | 6.3e-8 / 2.5e-6 |
| qC | 4 levels, every step | 1100 | 58 s | 2.7e-7 / 3.5e-6 |
| vA | 4 levels, every 2 | 1050 | 45 s | 3.0e-7 / 2.9e-6 |
| **vB** | **4 levels, every 4 (default of the script)** | **1150** | **45.5 s** (22.6 s on 7 ranks) | 1.8e-7 / 2.9e-6 |
| wA | 4 levels, every 8 | 1400 | 54 s | 4.7e-7 / 4.7e-6 |
| wB | 4 levels, every 16 | diverged / stalled (res_H 7e-3) | - | - |
| wC | 3 levels, every 4 | 1400 | 56 s | 6.3e-7 / 5.0e-6 |
| wD | 4 levels, every 4, `mg_coarse=1` | 1250 | 50 s | 1.2e-7 / 2.8e-6 |
| vC | 4 levels, every 4, `mg_start=50` | 1550 | 91 s | 4.7e-7 / 1.9e-5 |
| qD | 5 levels | diverged (res_H 6e-2) | - | - |
| vD | 4 levels, `mg_pre=mg_post=2`, `mg_coarse=4` | NaN at step ≈ 300 (extrapolation too long) | - | - |

Reading: the multigrid correction halves the number of iterations (1150 vs ≈ 2300) and about halves the wall time at a cost of 0–30 % per step (every-4 cycle: ≈ 0 %, every-step: ≈ +28 %).
`res_H` is lowered by two orders at equal iteration count (2000 it: 5e-5 → 2e-7), **`res_B` improves much less** (9e-5 → 3e-6 at 2000 it; its decay rate per iteration is the same,
only the starting level is lower): B is the bottleneck, see 6.4.

### 6.2 Step-ratio study of the induction (local pseudo time, no multigrid, 3000 iterations from the initial state, 3 ranks)

| `lts_b_ratio` R | res_H | res_B | remark |
|---|---|---|---|
| 1 (B on the global step) | 1.4e-4 | 3.8e-2 | stable, slow |
| 4 | 6.5e-6 | 9.0e-4 | |
| **16** | **1.1e-6** | **5.7e-6** | converged state: mass_H 3424.5, mass_Na 3.587 |
| 32 (with 4-level MG, 2000 it, run qB) | 6.3e-8 | 2.5e-6 | best |
| 64, with or without step smoothing κ = 0.15 / 0.3 | 2.5e-3 / 1e-3 / 5e-3 | | stagnates at 1e-3 or grows |
| 128 (+ κ=0.3) | 1.3e-1 | | diverges |
| ∞ (no cap) | – | – | mass_H grows from 3.3e3 to 5.6e4 by step 2000 and the run blows up (every CFL 1, 0.5, 0.3; with and without Hall/resistivity) |

Capping *both* cell and edge steps (`lts_max_ratio` 32 / 64 / 128, runs sA–sC) gives the same picture (32 ≈ same as R = 32, 64 stagnates, 128 diverges): the limit is the
**ratio of pseudo steps over the domain** (≈ 50), not the CT detail alone.  Smoothing the step field (`lts_smooth` κ = 0.2 / 0.5 / 1.0 without cap, κ = 0.15 / 0.3 with cap 64 / 128) did not help.

The converged states for R = 4, 16, 32 differ by up to 0.5 % in the H+ mass (3419.3 / 3424.5 / 3408.8): the exact fixed point is only recovered for uniform steps (section 2), the state
differs from the true steady state by an error concentrated where the step varies.  R = 1 is exact but slow.  The steady state of the converged runs is independent of the number of ranks and of the hierarchy (section 5).

### 6.3 Refined mesh (≈ 125 k cells, 7 ranks, start from the coarse steady state of run qB mapped by `mpcns_pre task=refine`, 3000 iterations, `lts_b_ratio=32`)

| | LTS only (fA) | LTS + MG, 4 levels, every 4, start 200 (fB) |
|---|---|---|
| clusters per level | – | 120234 → 29851 → 7336 → 1786 |
| res_H < 3e-3 at | step 900 (136 s) | step 500 (89 s) |
| res_H < 1e-3 at | not reached | step 1200 (214 s) |
| res_H < 5e-4 at | not reached | step 1400 (250 s) |
| res_B < 2e-2 at | step 900 | step 1300 |
| at step 3000: res_H / res_Na / res_B | 1.5e-3 / 5.7e-2 / 1.3e-2 | 3.2e-4 / 1.9e-2 / 7.0e-3 |
| wall for 3000 iterations | 454 s (151 ms/step) | 536 s (179 ms/step; +18 %) |

(Caveat: the starting coarse state was the R = 32 state of run qB, which is 15–70 % off the true steady state (6.5); the refined-mesh runs also use R = 32, so their fixed point is polluted in the same way — the comparison LTS vs MG is fair, the absolute state is not accurate.)

On the refined mesh the LTS iteration alone **stagnates** at res_H ≈ 2e-3 (the Na+ loading of the resolved exosphere, and a slowly relaxing near-planet region) while the multigrid
iteration keeps decreasing; neither has converged in 3000 iterations (Na+ mass is still moving: 3.1 → 5.5 → 5.0 for MG).  A refined-mesh run needs 1e4 or more iterations; the cost per iteration is
dominated by the GMRES of the semi-implicit induction (15–40 iterations per step at the large edge steps, ILU is not available with per-edge steps).

### 6.4 What limits the convergence now
1. **B.** `res_B` falls ≈ ×2.7 per 300 iterations on the coarse mesh and ≈ ×1.5 per 1000 on the refined mesh with or without MG.  The edge step is capped at R·τ_min (R ≤ 32) while the fluid cells far from the planet could use 100–300 τ_min.
   B has no multigrid correction.  A proper B correction needs a coarse CT hierarchy (coarse face fluxes = sums of the fine face fluxes over the children preserve $d_2B=0$ exactly; the nested octree refinement of `mpcns_pre`
   is a natural hierarchy, see `transfer_B`/`Transfer.cpp`), or an implicit treatment of the ideal induction.
2. **Hierarchy depth.** 5 levels diverge, more smoothing diverges: the coarse operator (first-order, frozen B, uniformly distributed correction) is a poor model of the near-planet physics with its strongly localised sources;
   the safeguard and the "start after 200 iterations" are needed.
3. **GMRES cost** of the implicit induction with large steps (up to 40 iterations/step) is now the largest per-step cost on the refined mesh.

### 6.5 Time-accurate reference and the accuracy of the pseudo-time fixed point (important)
Global-step run of the same coarse case (`steady=0`, 12 ranks, `time_order=1`, `second_order=0`, 5.4 ms/step, 36750 steps = 199 s to t = 75; run `base75`, fields `base75f`):

| t | 8.5 | 16.5 | 24.6 | 32.7 | 40.9 | 49.0 | 57.2 | 65.3 | 73.5 |
|---|---|---|---|---|---|---|---|---|---|
| res_H | 3.4e-2 | 6.0e-3 | 1.8e-3 | 3.7e-4 | 5.2e-5 | 9.4e-6 | 2.3e-6 | 5.6e-7 | 1.4e-7 |
| res_B | 4.2e-2 | 1.3e-2 | 2.3e-3 | 4.1e-4 | 8.0e-5 | 1.5e-5 | 3.4e-6 | 9.8e-7 | 6.4e-7 |

It converges (an e-fold per ≈ 6.6 time units = 3300 steps) to the **true steady state of the discrete model**: mass_H = 3419.05, mass_Na = 3.6613, min rho_H = 0.0433;
tolerance (res_H < 1e-5, res_B < 1e-3 or 1e-4) at step 24000 (res_H is the binding one; logged every 500 steps).  (An earlier statement in this log that it "never converges" was wrong: it had only been run to t = 25.)

The pseudo-time (LTS) iteration converges much faster **but to a slightly different state when the edge step is larger than the global one** (section 2: fixed point $d_1(TE)=0$).  Error of the converged LTS+MG state against the
time-accurate steady state (relative L2 over the fluid cells; B = induced face flux), coarse mesh, 3 ranks, 4-level MG every 4 steps:

| `lts_b_ratio` R | iterations to res_H<1e-5, res_B<1e-3 | rho_H | rho_Na | p_H | B flux | mass_H / mass_Na (true 3419.05 / 3.6613) |
|---|---|---|---|---|---|---|
| 1 | not reached (res_B 3.9e-2 at 3000) | 1.6e-2 | 2.0e-2 | 3.8e-2 | 0.51 (unconverged) | 3427.3 / 3.723 |
| 2 | not reached (res_B 1.9e-2 at 3000) | 1.0e-2 | 1.0e-2 | 1.8e-2 | 0.17 (unconverged) | 3419.0 / 3.560 |
| **4** | **2900** | **1.6e-2** | **3.7e-2** | **3.7e-2** | **3.3e-2** | 3419.3 / 3.622 |
| 8 | 1600 | 5.8e-2 | 0.19 | 0.14 | 0.12 | 3421.5 / 3.777 |
| 16 (LTS only) | 2500 | 0.11 | 0.32 | 0.27 | 0.29 | 3424.5 / 3.587 |
| 32 | 1150 (to res_B < 1e-4) | 0.15 | 0.40 | 0.32 | 0.71 | 3408.7 / 3.127 |

So: **R = 32 (the fastest) converges to a state that is 15–70 % off**; R = 4 is the largest value with an error of a few percent; only R = 1 is exact but then B converges slowly.
Cost to a tolerance (res_H < 1e-5, res_B < 1e-3), core-seconds on this node: time-accurate 24000 steps × 65 ms = 1560 core-s;
LTS + MG, R = 4: 2900 iterations × 90 ms = 260 core-s (**≈ 8× fewer iterations, 6× cheaper**), R = 8: 1600 × 95 ms = 150 core-s (15× fewer iterations, 10× cheaper, but 6–19 % error), R = 32: 1150 × 60 ms = 70 core-s (21× / 22×, 15–70 % error).
Use R = 32 only to produce a fast first approximation (e.g. as the coarse solution for the refinement indicators), R = 4 for a few-percent steady state, a time-accurate run for the exact one.

### 6.6 Multigrid at R = 4, and annealing (coarse mesh, 3 ranks, 4000 iterations)
| run | configuration | iterations to (res_H<1e-5, res_B<1e-3) | final res_Na | error vs true steady (rho_H / rho_Na / B) |
|---|---|---|---|---|
| eN4 | LTS only, R = 4 | 3000 | 4.8e-6 | 1.6e-2 / 3.7e-2 / 3.0e-2 |
| eR4c | + MG 4 levels, every 4 | 2900 | **1.7e-2** (floor) | 1.6e-2 / 3.7e-2 / 3.0e-2 |
| eR4b | + MG every 2 | not reached in 4000 (res_H 2.2e-5) | 2.2e-2 | same |
| gA / gB / gC | + MG with `mg_relax` 0.5 / 0.25 / `mg_growth` 1.2 | 2900 | 1.7e-5 / 3.1e-6 / 2.6e-5 | same |
| eN2 | LTS only, R = 2 | not reached (res_B 6e-3) | 9e-5 | 1.0e-2 / 7.0e-3 / 5.9e-2 |
| hA / hB | restart from the R = 32 state with R = 4 (annealing), without / with MG | hA: not reached, hB: 3000 | 3e-5 / 2e-5 | 2.6e-2 / 4.1e-2 / **0.47** |

* At R = 4 the number of iterations is the same with and without multigrid (B-limited); with `mg_relax = 1` the Na+ residual stalls at 1.7e-2 (the coarse correction and the fine step alternate); `mg_relax ≤ 0.5` or `mg_growth` 1.2 removes the floor at no loss.
* Annealing (fast R = 32 first, then R = 4) does **not** pay: the 71 % B-flux error of the R = 32 state decays only at the B pseudo-time rate (still 47 % after 3000 iterations at R = 4).

### 6.7 Prototype: implicit ideal induction (`implicit_ideal=1`, negative result)
Idea: with u frozen the ideal EMF (central average + transverse upwind) is linear in the induced flux, so $L(\Phi)$ joins $F$ in $(I+F d_1T)E=\dots$ (GMRES, Jacobi), with a *uniform* edge step $T=R\,\tau_{\min}$ (exact fixed point $d_1E=0$).
Code: `ImplicitInduction::set_ideal_state / ideal_linear`, `Config::implicit_ideal`; default off.  Coarse mesh, 3 ranks, local fluid steps, no MG:

| R | GMRES iterations / step (Jacobi) | outcome |
|---|---|---|
| 4 | 36–47 (explicit R = 4 needs ~2) | 3000 iterations took 15 min (explicit: 100 s): `res_H` stalls at 3e-3, `res_B` at 0.2 (limit cycle); the state is 1.3 % / 1.1 % / 2.7 % / 5.7 % off (rho_H / rho_Na / p / B), no better than explicit R = 4 |
| 16 | about 115 (rtol 1e-4, 500 max) | far too expensive |
| 64 | more than 260, does not converge | – |
| 256 | 100 (cap), residual 1 | – |

Two reasons: (i) without a real preconditioner the GMRES cost grows roughly with R (advective operator); (ii) freezing u makes the implicit step a Picard iteration of the u–B (Alfvén) coupling, which does not converge for large steps.  My earlier recommendation of this route was therefore wrong in practice; a usable implicit scheme would have to couple u and B (or use a much stronger preconditioner) and is not a small change.

## 7. How to run
```
export PATH=/workspace/software/openmpi-5.0.11/bin:$PATH
# coarse steady solution (mesh of 02_coarse_run.sh), 12 ranks:           scripts/05_steady_run.sh coarse 12 5000     (MODE=fast: R=32 + 4-level MG, default;  MODE=accurate: R=4, no MG)
# refine from it (03_refine.sh), then the refined steady run:           scripts/05_steady_run.sh fine 56 20000     (RESTART=1 continues)
# by hand:
mpirun -np 12 MercuryUns setup.txt mercury.txt geometry_dir=... output_dir=... steady=1 second_order=0 lts_b_ratio=32 \
    mg_levels=4 mg_cycle_every=4 mg_start=200 steady_tol=1e-7 max_Nstep=5000 log_every=100
python3 tools/conv_report.py run.log        # convergence table
python3 tools/compare_states.py <geometry_dir> ref/fields_final.unsb run/DATA/fields_final.unsb   # relative L2 error against a reference state
```
`MG_LEVELS=1` (script) / `mg_levels=1` switches the multigrid off (local steps only); `steady=0` is the old time-accurate solver.  The log line gets `res_H res_Na res_B res_Bsolid dt_max gmres_it gmres_res`;
at the end of a multigrid run the hierarchy sizes are printed once (`multigrid: 4 levels, clusters per level: ...`).

## 8. Limitations and next steps
* **Accuracy of the fast mode.**  The per-edge step changes the fixed point (section 2, table 6.5): the fast mode (R = 32) is a first approximation only.  The obvious cure, an implicit ideal induction with a uniform B step, was prototyped and **did not work** (section 6.7); what remains open is a B hierarchy (nested octree) or an implicit scheme that also couples u.
* The multigrid corrects the fluid only; B is frozen on the coarse levels and converges at its own pseudo-time rate (the slowest quantity, 6.4).  Next: B hierarchy from the nested octree (restriction of face fluxes = sum over the children faces; the prolongation is `transfer_B`), or implicit ideal induction so that the edge step can follow the cell step.
* The step ratio of the pseudo time is limited to ≈ 32–50 over the domain (section 6.2); the reason (beyond the CT term of section 2) is not fully understood.  A time-accurate-consistent multi-rate scheme or a Rosenbrock/implicit treatment of the Alfvén coupling would remove the limit.
* First-order coarse operator and piecewise-constant prolongation (so the convergence factor per cycle is modest); linear prolongation or a Galerkin-type correction are options.
* Only first-order fine discretisation was used (`second_order=0`): with `second_order=1` the fine rate comes from the MUSCL operator and the coarse levels are first order, which is the usual FAS combination but it was not run.
* The converged state of the pseudo-time iteration is the steady state of the model; the physical t = 25 state of the time-accurate run is not steady (Na+ is still loading), so the two answers differ by construction.
* Agglomeration is algebraic (heavy-face matching); it ignores the direction of B (anisotropy of the fast waves).  A weight that includes the Alfvén anisotropy might allow deeper hierarchies.
* The 56-rank production mesh has not been run in steady mode (only 15 test processes were allowed).

## 9. Development log (chronological)
1. **Baseline.**  Coarse mesh 24192 cells, 12 ranks, time-accurate `time_order=1`, `second_order=0`: 12500 steps to t = 25, `res_H` 1.8e-3, `res_B` 2.3e-3; mass_H 3399.8, min rho_H 0.043 (identical to the numbers in `MercuryRun.md`).
2. **Residual monitoring** added to the log (`res_H`, `res_Na`, `res_B`); first version of `res_B` used $d_1E$.
3. **Local pseudo time, version 0:** `cell_time_steps` (per-cell CFL), per-cell dt in `add_sources`/`gyrate`, `faraday_update_edge_dt` (B −= d1(dt_e E)), `ImplicitInduction::solve(..., dt_edge)` ($I + F d_1 T$).
   Result: with the unlimited step (ratio ≈ 250) the mass of H+ grows from 3.3e3 to 5.6e4 within 2000 iterations (CFL 1, 0.5, 0.3; with and without Hall/resistivity; also from the converged global-step state).
4. **Checks that the implementation is right:** GMRES converges (1e-4) in all runs; `lts_max_ratio=1` reproduces the global-step run (mass equal to 6 digits, later `test_pseudo_time` B to 2e-16).
5. **Diagnosis:** capping only the edge step (`lts_b_ratio`) stabilises the run: R = 8 → `res_H` falls to 2e-6.  The B residual measured with $d_1E$ stays at 0.2 although B is converged: the debug print of the worst face showed the
   actual update $d_1(TE)$ is 300× smaller than $\tau\,d_1E$.  Conclusion (section 2): the fixed point of the LTS-CT iteration is $d_1(TE)=0$, so $d_1E$ is the wrong residual and an unlimited per-edge step changes the steady state.
   `res_B` redefined as $\Delta\Phi/\tau_f$, `res_Bsolid` added.
6. **Step-ratio table** (6.2): R = 1, 4, 16 converge, 64 stagnates, ∞ blows up.  Step-field smoothing `lts_smooth` tried (0.2–1.0 uncapped, 0.15–0.3 with cap 64/128): no gain.  Cap of both steps `lts_max_ratio`: limit ≈ 50.
7. **Multigrid v0** (agglomeration, FAS, B frozen, per-level pseudo step with exact gyration and forcing outside the splitting).  First tests from the initial state: 2 and 3 levels run, 4 levels give `invalid time step` at the first correction
   (the correction injects up to 30 fine steps of Na+ production at once); `mg_start` (start after the transient) and the growth cap `mg_growth` added; 3 levels with start 300: `res_H` 4e-7 at 1500 iterations vs 7e-5 without MG.  4 levels still diverged.
8. **Unit tests for the hierarchy** (`test_multigrid`).  Test B (level-0 flux = fine first-order RHS) failed by 35 %: boundary faces of one kind but different normals had been merged into one face (corner clusters).
   Fix: merge only faces with $\cos>0.9$.  After the fix 4 levels converge (1100 iterations to tolerance).
9. **Parameter study** (6.1): cycle frequency (every 4 steps is as good as every step and cheaper), start iteration (200 better than 50), smoothing counts (more smoothing diverges), 3 vs 4 vs 5 levels.
10. **Safety:** `step_pseudo` throws `invalid time step (pseudo-time iteration diverged)` if no finite step exists (a diverged run used to continue with NaN).
11. **Refined mesh** test (7 ranks, 125 k cells): 6.3.  np independence: 6.1/5.
12. **Accuracy check against the time-accurate steady state** (section 6.5): the global-step run, continued to t = 75, converges (res 1e-7) to mass_Na 3.6613; the LTS fixed points differ by 3 % (R = 4) to 15–70 % (R = 32); earlier statements that the global-step run "never converges" were wrong and are corrected; annealing R tested (6.6).  `tools/compare_states.py` added.
13. **Implicit ideal induction prototype** (6.7): works mechanically, but GMRES cost grows with R and the frozen-u iteration does not converge; kept as an opt-in experimental switch.
14. Files added/changed: `MercuryUns/multigrid/*`, `driver/Solver.*` (`step_pseudo`, `mg_correct_`, residuals), `fluid/FluidRHS.*` (`cell_time_steps`, `add_sources_n`), `physics/Gyration.*` (`gyrate_n`, per-cell steps), `induction/Faraday.*` (`faraday_update_edge_dt`),
    `induction/ImplicitInduction.*` (edge steps), `Config.*`, `tests/test_multigrid.cpp`, `tests/test_pseudo_time.cpp`, `tools/conv_report.py`, `scripts/05_steady_run.sh`, `PLAN.md` (M9).

15. **Two attempts that were tried and reverted (not in the code).** (a) Anderson acceleration (depth 5–20, weighted inner product, on the whole state, B step R = 1): worse than the plain iteration (res_H 1e-3 instead of 3e-5 after 3000 iterations; the pseudo-step map is non-smooth because of floors, limiters and the varying step field).
    (b) Projected B update (desired increment −τ_f d1E projected on range(d1) by a few Jacobi-CG iterations of the curl-curl system; exact fixed point and exact d2 = 0 for any iteration count) was written but abandoned before the first run, by decision: the next approach is a curl-curl multigrid for B on a nested hierarchy (section 8).

## 10. References
* A. Jameson, W. Schmidt, E. Turkel, *Numerical solution of the Euler equations by finite volume methods using Runge-Kutta time-stepping schemes*, AIAA 81-1259 (local time stepping, multigrid).
* A. Brandt, *Multi-level adaptive solutions to boundary-value problems*, Math. Comp. 31 (1977) 333 (FAS).
* W. A. Mulder, *A new multigrid approach to convection problems*, J. Comput. Phys. 83 (1989) 303 (directional coarsening for convection).
* D. J. Mavriplis, V. Venkatakrishnan, *A 3D agglomeration multigrid solver for the Reynolds-averaged Navier-Stokes equations on unstructured meshes*, Int. J. Numer. Methods Fluids 23 (1996) 527.
* C. R. Evans, J. F. Hawley, *Simulation of magnetohydrodynamic flows: a constrained transport method*, ApJ 332 (1988) 659 (CT).

## 11. B multigrid (curl-curl, nested hierarchy): preparation

Direction decided with the user (2026-11-05): no per-edge local step for B (it changes the fixed point), no projection, no Anderson.  B is relaxed on every level with the *same*
uniform physical step and the coarse levels correct the long waves; restriction of the face flux = sum of the child face fluxes (exact discrete Stokes, $d_2$ preserved), prolongation = `transfer_B`;
fixed point stays $d_1E=0$.  Hierarchy: route 1 = structured coarsening of the cubed-sphere base mesh (2:1 in `cs_ncube`, `cs_nr`, `cs_nsolid`, same node lattice) + the octree levels of `mpcns_pre task=refine`.
All descendants of a coarse cell have to sit on the rank that owns it, hence a partition that is aligned to the hierarchy.

### 11.1 Partition control (done): `partition_mg_levels = K`
`mpcns_pre ... partition_mg_levels=K cs_ncube=M cs_nr=N cs_nsolid=S` (generate and refine): the base cells are grouped into blocks of $2^K\times2^K\times2^K$ (needs $M,N,S$ divisible by $2^K$), all cells of a block
(for `refine`: all fine cells whose base cell is in the block) stay on one rank; the blocks are partitioned on the quotient graph (`partition_cells_grouped`; weights = cost model, Hilbert + greedy refinement as before).
K = 0 is the old partition.  The code checks that every block is connected (guards the assumed generator cell order); `tests/test_partition_grouped.cpp`.
Load imbalance (max load / mean load) and cut faces, Mercury mesh (M = 12, N = 28, S = 4 → K ≤ 2: blocks of 8 or 64 cells, 3024 / 378 blocks on the base mesh):

| mesh | ranks | K = 0 (normal) | K = 1 (8-cell blocks) | K = 2 (64-cell blocks) |
|---|---|---|---|---|
| base 24192 cells | 7 | 1.004 / 6120 | 1.028 / 5480 | 1.027 / 4864 |
| | 12 | 1.013 / 7255 | 1.028 / 6212 | 1.020 / 6048 |
| | 28 | 1.029 / 10055 | 1.028 / 9132 | 1.062 / 8976 |
| | 56 | 1.030 / 12558 | 1.030 / 11552 | **1.143** / 11120 |
| refined 129738 cells (levels 0–2, the production-like setup) | 7 | 1.004 / 17588 | 1.025 / 11927 | 1.069 / 12954 |
| | 12 | 1.007 / 22587 | 1.029 / 16819 | 1.073 / 20620 |
| | 28 | 1.016 / 31039 | 1.029 / 25097 | **1.184** / 27434 |
| | 56 | 1.030 / 40739 | **1.109** / 34368 | **1.606** / 33283 |
(entries: imbalance / cut faces)

Reading: aligning to 8-cell blocks (K = 1, one multigrid level below the base mesh) costs ≤ 3 % imbalance up to 28 ranks (11 % at 56 ranks, 2300 cells/rank) and **reduces** the number of cut faces by 15–30 %
(blocky domains); 64-cell blocks (K = 2) cost 7 % at 7–12 ranks, 18 % at 28 ranks and 61 % at 56 ranks on the refined mesh, because a refined 64-cell block carries up to 64 × 64 fine cells and the
granularity (≈ 350 blocks of very unequal weight) is too coarse.  Recommendation: K ≤ 1 for the production partition (56 ranks), K = 2 only up to ≈ 12–28 ranks; the granularity problem would be
reduced by a finer weight model of the blocks (not tried) or by letting the refined region leave the block-aligned rule (levels below the base mesh then live on fewer ranks).

### 11.2 Two-level B prototype with frozen fluid (`mg_b_proto`, run 2026-11-05)

Program: `MercuryUns/bproto.cpp` (serial, one rank), new `Solver::b_dt / b_smooth / b_rate` (B-only step; fluid frozen; the same induction step as the pseudo-time solver — ideal + ambipolar EMF, implicit
resistive shell — with **one uniform step per level**, so the fixed point is exactly $d_1E=0$).  Hall term off (it costs 2× and is a ~1 % effect here).
Levels: coarse $C$ = base mesh (3024 cells, `cs_ncube=6, cs_nr=14, cs_nsolid=2`, stretch 1.15), fine $F$ = uniform 2×2×2 octree refinement of $C$ with spherical snap (24192 cells).
The fluid state is the pseudo-time steady state of $C$ (2500 LTS iterations, res_H 1e-10), prolonged piecewise constant to $F$ and kept frozen on both levels; induced B starts from zero on both.
FAS cycle for the flux (section 4.3 transposed to $\Phi$): $\nu_1$ fine smoothing steps; fine rate $r_F=(S_F(\Phi)-\Phi)/\Delta\tau_F$; $\Phi_C^0=R\Phi_F$ with $R$ = **sum of the child face fluxes** (sign by the area vectors);
forcing $g=Rr_F-r_C(\Phi_C^0)$ (it is $d_2$-free because $R$ commutes with $d_2$); $\nu_c$ coarse steps $\Phi_C\leftarrow S_C(\Phi_C)+\Delta\tau_C g$ with the coarse CFL step ($\Delta\tau_C/\Delta\tau_F=2.68$);
correction $P(\Phi_C-\Phi_C^0)$ with $P$ = `transfer_B` (exact $d_2=0$ on every fine cell, $RP=\mathrm{id}$), added with relaxation $\omega$.  Residual `res_B` = $\|d_1E\|/\|\Phi_{tot}\|$ of the last fine step.

Result (fine steps = smoothing steps on $F$; wall time on one core, includes the 1 fine + 1 coarse rate evaluation and one `transfer_B` per cycle; res_B at the same number of fine steps):

| fine steps | plain (uniform step, R = 1) | mgA ($\nu_1=5,\nu_c=5$) | **mgB ($\nu_1=10,\nu_c=10$)** | mgC (5, 5, $\omega$=0.5) | mgD (5, 20) |
|---|---|---|---|---|---|
| 200 | 58.9 s / 0.1827 | 119.8 s / 0.1815 | 98.7 s / 0.1814 | 111.6 s / 0.1767 | 139.9 s / 0.1950 |
| 600 | 157.6 s / 0.1610 | 331.5 s / 0.1475 | 273.5 s / 0.1478 | 312.4 s / 0.1568 | diverges (6.0) |
| 1000 | 255.1 s / 0.1443 | 512.9 s / 0.1172 | 427.7 s / 0.1167 | 490.6 s / 0.1234 | |
| 1400 | 349.2 s / 0.1257 | 652.8 s / 0.1035 | 563.8 s / 0.1029 | 625.1 s / 0.1059 | |
| 2000 | 481.4 s / 0.1055 | | 725.7 s / 0.0894 | | |

* Same residual after fewer fine steps: mgB reaches 0.1055 after ≈ 1400 fine steps instead of 2000 (**≈ 1.4–1.7× fewer iterations**, 1.9× at 800–1000 steps), but the wall time is **not** reduced (≈ 540 s vs 481 s):
  the two rate evaluations (each a full step) and `transfer_B` cost as much as the saved steps.  Coarse relaxation 0.5 is worse, more coarse steps ($\nu_c$ = 20) diverge.
* The iteration is still far from converged after 2000 steps (res_B 0.1 → needs ~1e-6): with a uniform step the fixed-point iteration for B in a frozen flow crawls.

**Why the gain is bounded (analysis).**  With the fluid frozen the B equation is hyperbolic (advection by $u$ with the stiff resistive shell treated implicitly): information moves one cell per step at most, so the number of
steps to renew the field across the domain is $L/(u\,\Delta\tau)$ with $\Delta\tau=\Delta\tau_{\min}$ set by the smallest, most strongly magnetised cells (here 2.2e-3: ≈ 6000 steps for L = 13, u = 1).
A coarse level lets the long waves travel with its own CFL step, so the gain per level is at most $\Delta\tau_C/\Delta\tau_F$ (= 2.68 here; measured 1.4–1.9×).  This is the known limit of multigrid for
hyperbolic problems; it is the same reason why the local pseudo time was so effective for the *fluid* (steps ×250 in the far field) and why the per-edge step is what made B fast but wrong.
A hierarchy of $m$ levels can give at most $\sim 2.7^{\,m-1}$ (≈ 20× for 4 levels), but each extra level needs the structured coarsening of the base mesh and costs the same overheads.

**Conclusion of the two-level experiment:** the uniform-step B multigrid is exact (fixed point $d_1E=0$, $d_2=0$ preserved to round-off through $R$ and `transfer_B`) but gives ≈ 1.5× fewer iterations and no wall-time gain with two levels;
the expected benefit with more levels is bounded by the CFL ratio of the levels.  Cheaper variants to try next: use the rate of the last smoothing step instead of a separate evaluation, a cheaper $P$ (local, linear) instead of the CG inside `transfer_B`, 3–4 levels.

### 11.3 Multigrid on / off with the Hall term ON (matched runs, 2026-11-05)
Decision of the user: the Hall term must stay on.  All steady-mode runs of sections 6.1–6.6 used `mercury.txt` (Hall on); only the B prototype of 11.2 had it off.  Matched comparison, coarse mesh, 3 ranks, Hall on,
`mg_levels=4`, cycle every 4 steps, start 200 (tolerance res_H < 1e-5, res_B < 1e-3):

| mode | multigrid | iterations | wall | final res_H / res_B (3000 it) |
|---|---|---|---|---|
| fast (R = 32) | off | 2700 | 90 s | 4.1e-6 / 1.3e-5 |
| fast (R = 32) | **on** | **1100** | **38 s** | 1.9e-7 / 7.2e-7 |
| accurate (R = 4) | off | 2950 | 85 s | 6.5e-6 / 9.0e-4 |
| accurate (R = 4) | on (`mg_relax=0.5`) | 2900 | 94 s | 1.2e-6 / 7.7e-4 |

Refined mesh (129 k cells, 7 ranks, R = 32, 3000 iterations): res_H 1.5e-3 → 3.2e-4, wall +18 % (6.3).
Two-level B multigrid (fluid frozen), Hall ON, 800 fine steps: plain 214 s / res_B 0.1257, two-level (ν1=10, νc=10) 326 s / 0.1069; plain needs ≈ 1100 steps (≈ 300 s) for 0.107: 1.4× fewer iterations, 1.1× more wall time.
Summary: the fluid multigrid is worth it in the fast mode (2.4× fewer iterations, 2.3× less wall time; more on refined meshes), useless in the accurate mode (B-limited); the B multigrid gives no wall-time gain.
