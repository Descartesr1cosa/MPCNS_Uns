# Ideal EMF on edges of arbitrary valence

Code: `Common/induction/IdealEMF.{h,cpp}`, weights `Geometry::ef_weight`, test `tests/test_ideal_emf.cpp`.

## Scheme
Cartesian form of the legacy Mercury scheme (`MercurySolver_MagExplicit_Ideal.cpp`, Gardiner–Stone type UCT).

**Face candidate** (n = canonical face normal pointing from cell L to cell R, L = cell whose `c2f` sign is +1):
`F_X = (u_X·n) B_X^tot − u_X B_f^tot,n`,  `F = ½(F_L+F_R) − ½ a_f (B_R^ind − B_L^ind)`,  `a_f = max(|u_L·n|,|u_R·n|)`,
`E_f = −n × F` (tangential to the face).  A boundary face uses its single cell, without dissipation.

**Edge value** (line integral along the stored edge orientation `t_e`):

    E_e = Σ_f w_{e,f} E_f·t_e  +  2 C L_e Σ_f ω_f a_{e,f} (t̂_e·(m_f × n_f)) B_n,f

- `w_{e,f}`: share of face f in the dual face of e (pieces of both cells adjacent to f), Σ_f w = 1; valence 4 → 1/4.
- `m_f`: unit in-face direction pointing away from the edge; `a_{e,f} = max_{c∋e} |u_c·m_f|`; `B_n,f = Φ_f/|A_f|`.
- `ω_f = w_f (1 − q_f·λ)`, `q_f = a_f (t̂·(m_f×n_f)) n_f` projected on the plane normal to the edge, with λ the solution of the
  2×2 system `Σ_f w_f q_f q_fᵀ λ = Σ_f w_f q_f`.  This is the minimum change of the weights (in the norm Σ(Δω)²/w) for which
  `Σ_f ω_f q_f = 0`, i.e. **a uniform B produces no upwind term on any edge** (valence 3, 5, 6, 8, hanging nodes).
- In the legacy scheme the upwind term couples the two *opposite* faces across the edge; that pairing does not exist for
  other valences.  Without the correction a uniform B would get an O(a|B|) spurious term on irregular edges.
- On a regular hexahedral edge `g = Σ w q = 0` ⇒ `ω = w = 1/4` and the legacy formula is reproduced exactly.
- The 2×2 solve is done in an orthonormal basis of the plane normal to the edge.  (A regularised 3×3 Cramer solve loses all
  digits along the edge direction and made the result depend on the summation order — found and fixed while testing np-independence.)

## Evaluation and parallel use
Only the *owned* edges of a rank are evaluated: every cell and face around such an edge is local (vertex-neighbour ghost layer),
then `EdgeEMF::finish` copies the value to the alias edges.  Per step: `forward(Face, B)` → cell reconstruction of B → EMF →
`finish` → `faraday_update` on the faces of owned cells.  Cell fields (u, ...) must have been forwarded to ghost cells.
Other EMF terms (Hall, resistive, artificial) are added to the same edge array before `finish`; cell-based terms use
`EdgeEMF::add` from owned cells and the reduction path of `finish`.

## Verified (tests/test_ideal_emf.cpp, np = 1, 4, 7)
| property | result |
|---|---|
| regular hex, random u and B vs a transcription of the legacy formulas | 450 interior edges, max difference 4.0e-15 |
| uniform B, random u: upwind part | 5e-17 (box), 2e-14 (hanging), 3e-15…5e-17 for valence 3, 5, 6, 8 |
| uniform u, B on planar meshes: E = −(u×B)·dl | 1e-16 … 4e-15 |
| 20 CT steps, np = 1/4/7 vs serial | ≤ 9e-16 difference in B, max relative ∇·B ≤ 2e-15 |

## Known limits
- **Warped faces** (cubed-sphere spherical faces): the face candidate only determines the in-face components of E; the part
  `(n·E)(n·t)` of the edge projection is missing and the uniform-field error is first order: mean relative error of
  `−u×B` 2.8e-3 → 1.3e-3 → 5.9e-4 for ncube = 4, 8, 16 (the legacy scheme has the same property on curvilinear grids).
- **Boundary edges** (open face fan): the correction cannot satisfy the uniform-B constraint with fewer than three independent
  faces, so the upwind term is not exactly zero there (≈1e-12 on planar meshes, ≈1e-2 on the shell).  Boundary conditions
  will overwrite these edges; the legacy code relied on ghost cells instead.
- ω_f can become small or negative on strongly irregular edges; no limiter yet.  Stability on such edges is untested in a time loop with real dynamics.
- Background-field split: `phi_tot`/`B_cell_tot` vs `phi_ind`/`B_cell_ind` are separate inputs (legacy `Badd`); the solver does not provide Badd yet.
