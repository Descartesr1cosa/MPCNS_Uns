# Hall term — IMPLEMENTED (first order, semi-implicit); this file keeps the original design notes

Status 2026-11: the Hall EMF is implemented in `Common/induction/ImplicitInduction.*` (see "Semi-implicit induction" in `docs/MercuryRun.md`);
the stub `hall/HallEMF.*` is obsolete. Still open from the list below: the explicit whistler diffusion (not needed with the implicit treatment),
second-order time integration of the Hall term, sub-cycling (not needed).

# Hall term — deliberately left empty for now (`Common/induction/HallEMF.*`)

`add_hall_emf()` is a no-op and `hall_time_step_limit()` returns +inf.  The induction path
(`EdgeEMF` → `Exchanger::assemble_edges` → `faraday_update`) is complete and tested without it.

## Agreed design to implement later
1. **Form** (unchanged from Mercury): `E_hall = α (J×B)`, `α = hall_coef / n_e` with the existing
   `Hall_Num_Limiter`, radial taper `HallRadialTaper_` (tanh between `hall_taper_r_min/max`).
2. **Edge weights**: cell values `α_c J_c×B_c` are brought to an edge with the dual-face shares `w_{e,c}`
   (`Geometry::ec_weight`, Σ_c w = 1), projected on the edge: `E_e = Σ_c w_{e,c} (α J×B)_c · t_e`.
   Valence-4 edges give the legacy 1/4 average.  The term is added with `EdgeEMF::add` so it uses the same
   edge assembly and the same owner as the ideal EMF.
3. **J on edges**: `J_e = Σ_{f∋e} d1[f,e] (⋆2 B)_f` (transposed incidence, any valence); cell `J_c` by the
   same reconstruction used for `B_c`.
4. **Whistler-type diffusion**: legacy form `Cwh · a_w · jump(B_t)` with `a_w ∝ α|B|/h`, assembled with the
   same weights `w` so it applies to any edge valence.
5. **Time step**: `ω_max ≈ α|B|π²/h²` ⇒ `Δt_Hall = C_h h²/(α|B|)`, with `h_c = 6V/Σ|A|` and the shortest edge.
   Reported next to the CFL step.  Each refinement level halves the CFL step but divides the Hall step by 4.
6. **Relief** (CT is preserved because E stays an edge 1-cochain):
   A. B sub-cycling with fluid frozen; accumulate the time-integrated edge EMF and apply
      `B -= d1(Σ dt_s E_s)` once.
   B. only if profiling demands: region-local sub-cycling, or a semi-implicit Hall step
      (implicit midpoint + Newton–Krylov on edge dofs, PETSc).

Legacy references: `legacy/MPCNS_Mercury/Z4_Mercury/src/physics/MercurySolver_MagExplicit_Hall.cpp`,
`.../MercurySolver_Current.cpp`, `.../MercurySolver_Timestep.cpp`.
