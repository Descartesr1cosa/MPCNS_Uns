# Numerical methods of MPCNS_Uns — mathematical summary and references

Scope: what the code does, written as formulas, with the place in the code and the literature the technique comes from.
Everything is nondimensional (reference scales `U_ref, L_ref, B_ref, n_ref, T_ref`, see `MercuryUns/physics/Scales.cpp`).
Reference numbers `[R#]` are listed in section 11. Statements marked *(empirical)* were measured here, not proven.

Notation: cells $c$, faces $f$, edges $e$, nodes $n$; $\Phi_f$ = magnetic flux through face $f$ (canonical orientation); $E_e$ = line integral of the electric field along edge $e$;
$D=d_1$ the face–edge incidence matrix, $\Phi^{n+1}=\Phi^n-\Delta t\,DE$; $\star_k$ diagonal Hodge stars; $t_e$ edge vector, $A_f$ face area vector, $V_c$ volume.

---

## 1. Signed cell complex and discrete exterior calculus (`include/uns/complex`, `docs/Conventions.md`)

The mesh is stored only as signed incidence tables
$$e2n_{e,n}\in\{-1,+1\},\qquad f2e_{f,e}\in\{-1,0,+1\},\qquad c2f_{c,f}\in\{-1,0,+1\}.$$
The discrete operators are $d_0=e2n^T$ (nodes→edges), $d_1=f2e$ (edges→faces), $d_2=c2f$ (faces→cells), and
$$d_1d_0=0,\qquad d_2d_1=0 \tag{1.1}$$
hold **as integer identities of the tables** (checked by exact integer products, `count_nonzero_products`), independently of cell type, hanging nodes, partitioning and of the coordinates [R1,R2,R3].
Canonical orientation is derived from global node ids (edge: low→high id; face loop: starts at the smallest node, second node smaller than the last), so every rank agrees without communication; local copies carry an `alias_sign` relative to the owner.

Because $B$ is only ever changed by $\Phi\leftarrow\Phi-\Delta t\,DE$, with $E$ an arbitrary edge array,
$$(d_2\Phi^{n+1})_c=(d_2\Phi^n)_c-\Delta t\,(d_2d_1E)_c=(d_2\Phi^n)_c \tag{1.2}$$
exactly (up to floating-point round-off of the sums, measured $10^{-15}$–$10^{-12}$ relative). This is constrained transport [R4,R5] written on a general complex.

## 2. Metric (Hodge stars, dual weights) (`include/uns/geometry`)

Face centre $c_f$ = area-weighted centroid of the triangle fan of the loop; $A_f=\tfrac12\sum_i (v_i-c_f)\times(v_{i+1}-c_f)$ (well defined for warped faces).
Cell volume with cell centre $x_c$: $V_c=\tfrac13\sum_{f}\sum_{i} a_{f,i}\cdot(c_f-x_c)$ over the face-centre triangulation (watertight, independent of $x_c$).
Dual face of edge $e$ in cell $c$: triangles $(m_e,c_f,x_c)$ over the faces $f\ni e$ of $c$, with vector areas $\pi_{e,c,f}=\tfrac12 s\,(x_c-m_e)\times(c_f-m_e)$ oriented along $t_e$ ($s=\pm1$ from the incidence signs). Weights
$$w_{e,c}=\frac{\sum_{f}|\pi_{e,c,f}|}{\sum_{c',f}|\pi_{e,c',f}|},\qquad w_{e,f}=\frac{\sum_{c\ni f}|\pi_{e,c,f}|}{\sum_{c',f'}|\pi_{e,c',f'}|},\qquad \sum_c w_{e,c}=\sum_f w_{e,f}=1, \tag{2.1}$$
valid for any edge valence (valence 4 on a regular hexahedral grid gives $1/4$). Diagonal Hodge stars (all positive for a valid mesh) [R2,R6]:
$$\star_1(e)=\frac{|A^*_e|}{|t_e|},\quad \star_2(f)=\frac{\sum_{c\ni f}|c_f-x_c|}{|A_f|},\quad \star_3(c)=\frac1{V_c},\quad \star_0(n)=\text{dual volume}.$$
Current through the dual face of $e$: $J_e=(D^T\star_2\Phi)_e$ (circulation of $H=\star_2\Phi$ around the dual face).

## 3. Non-conforming refinement and divergence-free transfer (`src/uns/refine`)

**Mesh.** Every coarse hexahedron carries an octree. A lattice point $(a,b,c)\in[0,2^L]^3$ of a coarse cell is keyed by the coarse entity it lies on (node / edge / face / cell) and its integer coordinates in that entity's canonical frame, so shared points need no floating-point matching [R7]. Each leaf side is a dyadic square; the faces of the fine complex are the *common refinement* of the squares of the two sides of a plane: a face is a square that contains no smaller square of the opposite side, and the coarser neighbour lists the smaller squares as its own faces (a polyhedron). Every edge is split at all lattice nodes on it (found by recursive bisection through the node hash). Face balance (level difference ≤1 across faces) is enforced by iterating. Hanging nodes are thus ordinary vertices of the coarse neighbour and (1.1) holds by construction.

**B transfer.** Coarse flux $\Phi_F$ of coarse face $F$ is distributed on its fine tiles $T_j$:
$$\Phi_{T_j}=\hat\Phi_{T_j}+\frac{|A_{T_j}|}{\sum_k|A_{T_k}|}\Big(\Phi_F-\sum_k\hat\Phi_{T_k}\Big),\qquad \hat\Phi_T=\bar b\cdot A_T,$$
with $\bar b$ the mean of the cell vectors $b_c=\big(\sum_f A_f n_fn_f^T\big)^{-1}\sum_f A_f n_f\,\Phi^{out}_f/A_f$ (area-weighted least squares) of the adjacent coarse cells, so $\sum_j\Phi_{T_j}=\Phi_F$ exactly. Inside each refined coarse cell the interior fluxes start from $b_c\cdot A_f$ and receive the minimum-norm correction $\delta$ that makes every fine cell divergence free:
$$\min\sum_f\frac{\delta_f^2}{|A_f|}\ \ \text{s.t.}\ \ D_c\delta=-r,\ \ r=D_c\Phi^{0}\ \Rightarrow\ (D_cWD_c^T)\lambda=-r,\ \ \delta=WD_c^T\lambda,\ W=\operatorname{diag}|A_f|,$$
solved by Jacobi-preconditioned CG on the cell-graph Laplacian followed by a spanning-tree sweep that removes the remaining round-off residual (measured $|d_2B|/\sum|B_f|\le10^{-11}$). Cell fields are prolonged piecewise-constant by default (positive and conservative); a limited linear prolongation is optional.

## 4. Partitioning and load balance (`src/uns/partition`, `src/uns/parallel`, `MercuryUns/driver/Solver.cpp::write_load_profile`)

Graph: vertices = cells with weight $w_v=c_{\rm type}\,(1+\kappa\max(0,n_{\rm faces}-6))$ ($c_{\rm fluid}=1$, $c_{\rm solid}$ from the profile), edges = shared faces with weight $c_{\rm comm}$.
Partitioners: space-filling-curve ordering (Hilbert [R8], Morton) with a weighted prefix cut, $p(c)=\lfloor (W_{<c}+\tfrac12 w_c)/W\cdot N_p\rfloor$; weighted recursive coordinate bisection; greedy boundary refinement (move a boundary cell to the neighbour part with the largest cut gain if both loads stay within $(1\pm\text{tol})\bar w$) [R9 for the graph model].
Ownership: owner of an entity = lowest rank among the ranks of its incident cells; a rank stores its owned cells plus one vertex-neighbour layer, so every edge of an owned cell sees all its incident cells and faces locally.
Edge assembly (all EMF terms): partial sums over owned cells/faces are reduced to the owner and broadcast back ($E_e=\sum_{\rm ranks}E_e^{(r)}$), which also makes all copies of an edge bitwise equal.

Calibration from a sparse run. With per-rank compute time per step $\tau_r$ = (wall − communication/wait)/steps and owned fluid/solid cell counts $n^f_r,n^s_r$, the least-squares fit
$$\min_{a,b\ge0}\sum_r\big(\tau_r-a\,n^f_r-b\,n^s_r\big)^2,\qquad c_{\rm solid}=b/a,$$
and $c_{\rm comm}=\operatorname{median}_r\big(\text{wait}_r/n^{\rm ghost}_r\big)/a$ are written to `load_profile.json/.txt` and read by `mpcns_pre cost_profile=`. *(empirical: rank-to-rank differences are small, $r^2$ is low; it is a guide, not a fit to trust to 10 %.)*

## 5. Finite-volume fluid (`MercuryUns/fluid`)

Per species $s\in\{H^+,Na^+\}$, $U=(\rho,\rho u,E)$, $p=(\gamma-1)(E-\tfrac12\rho|u|^2)$:
$$\frac{dU_c}{dt}=-\frac1{V_c}\sum_f c2f_{cf}\,\mathcal F_f|A_f|+S_c,\qquad \mathcal F_f=\tfrac12(F(W_L)+F(W_R))\cdot n-\tfrac12\lambda_f(U_R-U_L),\ \ \lambda_f=\max_{L,R}(|u\cdot n|+c_s)$$
(Rusanov [R10]). Second order (optional; the production runs so far use first order): linear reconstruction of the primitives $W=(\rho,u,p)$ with the inverse-distance weighted least-squares gradient [R11]
$$\nabla W_c=\Big(\sum_k \omega_k d_kd_k^T\Big)^{-1}\sum_k\omega_k d_k(W_k-W_c),\quad \omega_k=|d_k|^{-2},$$
Barth–Jespersen limiter $\psi_c=\min_f\phi_f$ [R12] ($\phi_f$ from the neighbour min/max), face state $W_f=W_c+\psi_c\nabla W_c\cdot(x_f-x_c)$ with fall-back to $W_c$ if $\rho$ or $p$ would become negative. Boundary faces contribute a mirrored virtual neighbour so that boundary cells keep a full-rank stencil.
Boundary states: wall — $H^+$ reflected ($u\to u-2(u\cdot n)n$), $Na^+$ absorbing (ghost $\rho=p=10^{-12}$); inflow — free stream; outflow — zero gradient.

## 6. Two-ion model and source terms (`MercuryUns/physics`, `fluid/FluidRHS.cpp`)

Number densities $n_s=\rho_sM_{ref}/M_s$, $n_e=n_H+n_{Na}$, $n_{e,\rm eff}=\sqrt{n_e^2+n_{\rm floor}^2}$, $\chi_s=n_s/n_e$, $w_s=n_s/n_{e,\rm eff}$, $u^+=\chi_Hu_H+\chi_{Na}u_{Na}$.
Momentum sources (`inver_MA2`$=B_{ref}^2/(\mu_0\rho_{ref}U^2)$, `induce`$=qL B_{ref}n_{ref}/(\rho_{ref}U)$):
$$S^{\rm mom}_s=\underbrace{\text{induce}\,n_s(u_s-u^+)\times B}_{\text{gyration, §8.2}}+\text{inver\_MA2}\,w_s\,(J\times B)-w_s\nabla p_e-a_4\rho_su_sv_{st},\qquad p_e=p_H+p_{Na},\ J=\nabla\times B_{\rm ind},$$
energy source $u_s\cdot S^{\rm mom}$ (+ midpoint correction §8.3), Na$^+$ production $a_1q$ and internal energy $a_6qT_{n0}$ with $q=\nu\,n_{\rm Na,neutral}$ from the sodium exosphere model [R13] and a day/shadow photo-ionisation frequency $\nu\in\{s_{k1},s_{k2}\}$; temperature relaxation $a_5n_sb_1+a_4\rho_su_s^2b_2$ with $b_1=v_{st}(T_{n0}-T_s)m_s$, $b_2=m'v_{st}$. The model and the coefficient definitions follow the legacy MPCNS_Mercury source routine; multi-fluid/Hall-MHD references: [R14,R15].
Stiffness check: $a_4v_{st}\sim3\cdot10^{-4}$, relaxation $\sim10^{-4}$, production a bounded source, all $\ll1/\Delta t\sim10^4$ — treated explicitly.
Magnetic field: $B=B_{\rm ind}+B_{\rm add}$, $B_{\rm add}=B_{\rm IMF}+B_{\rm dipole}$; the dipole enters as the exactly divergence-free flux $\Phi_{\rm dip}=d_1(A\cdot dl)$ of the vector potential $A=m\times(r-r_0)/(|r-r_0|^2+a^2)^{3/2}$.

## 7. Electric field terms as edge 1-cochains

**7.1 Ideal EMF on edges of any valence** (`induction/IdealEMF.cpp`, `docs/IdealEMF.md`), a Cartesian generalisation of the Gardiner–Stone UCT scheme [R5,R16] used in the legacy code. Face candidate (normal $n$ from cell $L$ to $R$):
$$\mathcal G_f=\tfrac12(F_L+F_R)-\tfrac12a_f(B^{\rm ind}_R-B^{\rm ind}_L),\quad F_X=(u_X\!\cdot n)B_X-u_XB_{f,n},\quad E_f^{\rm face}=-n\times\mathcal G_f,$$
$$E_e=\sum_{f\ni e}w_{e,f}\,E_f^{\rm face}\!\cdot t_e\;+\;2C\,|t_e|\sum_{f\ni e}\omega_f\,a_{e,f}\,\sigma_f\,\frac{\Phi^{\rm ind}_f}{|A_f|},\qquad \sigma_f=\hat t_e\cdot(m_f\times n_f),$$
$m_f$ = unit in-face direction pointing away from the edge, $a_{e,f}=\max_{c\ni e}|u_c\cdot m_f|$, $C=0.5$. On a regular hexahedral edge the transverse (upwind) part couples the two opposite faces; for other valences the weights $w_f$ are corrected by the minimum change that makes a uniform field produce no upwind term:
$$\min\sum_f\frac{(\omega_f-w_f)^2}{w_f}\ \text{s.t.}\ \sum_f\omega_fq_f=0,\ \ q_f=a_{e,f}\sigma_fn_f^{\perp}\ \Rightarrow\ \omega_f=w_f(1-q_f\!\cdot\!\lambda),\ \ \Big(\sum_fw_fq_fq_f^T\Big)\lambda=\sum_fw_fq_f ,$$
solved as a $2\times2$ system in an orthonormal basis of the plane normal to $t_e$ (a regularised $3\times3$ Cramer solve loses all digits along $t_e$). Reduces exactly to the legacy formula on regular hexahedra (max difference $4\times10^{-15}$ on 450 edges).

**7.2 Ambipolar field** (`physics/Ambipolar.cpp`): $E_e^{\rm amb}=-p_{\rm coef}\,\dfrac{\tilde p_e(n_2)-\tilde p_e(n_1)}{n_{e,\rm eff}(e)}$, $\tilde p_e$ = volume-weighted nodal electron pressure, $n_{e,\rm eff}(e)=\sum_cw_{e,c}n_{e,\rm eff,c}$, $p_{\rm coef}=\rho_{ref}U/(qLB_{ref}n_{ref})$ — the exact $d_0$ of a nodal potential divided by an edge density.

**7.3 Resistivity** (Mercury interior), $\eta(r)=\eta_{\max}/(\mu_0UL)\cdot\tfrac12[\tanh\frac{r-r_1}{w}-\tanh\frac{r-r_2}{w}]$, $E^{\rm res}_e=\dfrac{\eta_e}{\star_1(e)}J_e=(M D^T\star_2\Phi)_e$, $M=\operatorname{diag}(\eta_e/\star_1(e))$ (current density $\times$ length: $J_e/|A^*_e|\cdot\eta\cdot|t_e|$). No artificial resistivity is used.

**7.4 Hall field**, $E^{H}=\dfrac{h_{\rm coef}}{n_{e,\rm eff}}\,\tau(r)\,J\times B$, $h_{\rm coef}=\dfrac{B_{ref}M_{ref}}{UqL\mu_0\rho_{ref}N_A}$, taper $\tau(r)=\tfrac12[1+\tanh(2(r-r_m)/\Delta r)]$ ($r_m$, $\Delta r$ from $r_{\min}=1.0$, $r_{\max}=1.012$); on edges: $E^H_e=\sum_cw_{e,c}\,\alpha_c\,(J_c\times B_c)\cdot t_e$, $\alpha_c=\tau_ch_{\rm coef}/n_{e,\rm eff,c}$ [R15,R17]. Whistler dispersion for $B_0\hat x$: $\omega=\alpha B_0k^2$.

## 8. Time integration

**8.1 `time_order=2`**: SSP-RK2 [R18] for $U$ and $\Phi$: $\Phi^{n+1}=\Phi^n-\tfrac{\Delta t}2D(E^n+E^{(1)})$, the average of two CT updates, so (1.2) is kept; resistivity as an implicit split step (§9 with $H=0$).

**8.2 Exact integration of the gyration term (Strang splitting [R19]).** The induce force is a rotation of the relative velocity $w=u_H-u_{Na}$ about $B$:
$$\dot w=\Omega_{\rm rel}\,w\times\hat b,\quad \Omega_{\rm rel}=(a_H\chi_{Na}+a_{Na}\chi_H)|B|,\ a_s=\text{induce}\,M_{ref}/M_s,$$
$$w(t)=w_\parallel+\cos\theta\,w_\perp+\sin\theta\,(w_\perp\times\hat b),\ \theta=\Omega_{\rm rel}t;\qquad \Delta u_H=f_H\Delta w,\ \Delta u_{Na}=-f_{Na}\Delta w,\ \ f_H=\frac{a_H\chi_{Na}}{a_H\chi_{Na}+a_{Na}\chi_H}=\frac{\rho_{Na}}{\rho_H+\rho_{Na}}.$$
Momentum and total kinetic energy are conserved; the energy is updated by $\tfrac12\rho(|u'|^2-|u|^2)$ so the thermal energy is untouched. Step: rotation $\Delta t/2$, explicit step, rotation $\Delta t/2$ (second order in the splitting).

**8.3 Work of the force with the midpoint velocity.** For an explicit step with momentum increment $\Delta u=\Delta t\,S/\rho$ the kinetic energy changes by $\rho u\cdot\Delta u+\tfrac12\rho|\Delta u|^2$; using $S_E=u\cdot S+\dfrac{\Delta t}{2\rho}|S|^2=u_{\rm mid}\cdot S$ makes the work equal to the kinetic-energy change, so the pressure $p=(\gamma-1)(E-\tfrac12\rho u^2)$ is not drained in low-β cells (β $\sim10^{-6}$ near the poles).

**8.4 Time step.** $\Delta t=\min_{c,s}\ \mathrm{CFL}\,V_c\Big/\sum_{f\subset c}(|u_s\cdot n_f|+c_{f,s})|A_f|$ with fast speed $c_f^2=\tfrac12\big[(c_s^2+v_A^2)+\sqrt{(c_s^2+v_A^2)^2-4c_s^2v_{A,n}^2}\big]$, $v_A^2=\text{inver\_MA2}\,|B|^2/\rho_{\rm tot}$ with the **total** density (floored at $10^{-2}$) because the Lorentz force acts on each species through its share $w_s$. Explicit whistler estimate (information only): $\Delta t_H=\min_c\, \ell_c^2\,n_{e,\rm eff}/(\pi^2\,h_{\rm coef}\,|B|)$ with $\ell_c=6V_c/\sum_f|A_f|$ and $h_{\rm coef}$ the Hall coefficient of §7.4 [R15].

## 9. Semi-implicit induction with $E$ as the unknown (`induction/ImplicitInduction.*`) — first order in time

**Formulation.** Split $E=E_{\rm expl}+F(\Phi^{n+1})$ with $E_{\rm expl}$ = ideal + ambipolar EMF at $t^n$ (and boundary values), and the implicit linear operator
$$F(\Phi)=M\,D^T\star_2\Phi+H(\Phi),\qquad H(\Phi)_e=\sum_{c\ni e}w_{e,c}\,\alpha_c\big(J_c[\Phi]\times B^0_c\big)\cdot t_e ,$$
$B^0$ = total cell field at $t^n$ (**linearised Hall**: $J^{n+1}\times B^n$) [R15,R20]; $J_c[\Phi]=\dfrac1{V_c}\sum_f s_{cf}A_f\times\tfrac12(b_c+b_{c'})$, $b_c=\sum_f R_{cf}\Phi_f$ the least-squares cell field of §3 (precomputed matrix $R$). Substituting $\Phi^{n+1}=\Phi^n-\Delta tDE$:
$$\boxed{\ (I+\Delta t\,F\,D)\,E=E_{\rm expl}+F(\Phi^n)\ },\qquad \Phi^{n+1}=\Phi^n-\Delta t\,DE. \tag{9.1}$$
The unknown is the edge 1-cochain $E$, hence (1.2) holds **whatever the accuracy of the linear solve**; this is why $E$, not $B$, is the implicit variable. Imposed inflow EMFs are Dirichlet unknowns removed from the system (their rows are identity and the Krylov corrections vanish there); wall edges ($E_{\rm expl}=0$) stay unknowns so the tangential $\eta J$ appears. The sum over cells is evaluated by owned cells and reduced through the edge owners (one edge exchange per matrix–vector product, `assemble_edges` with two components).

**Why this method.** Backward Euler for the resistive part is unconditionally stable; with the frozen-$B^0$ Hall operator the system is linear (no Newton iteration, in contrast to a fully nonlinear implicit/JFNK treatment [R21]) and unconditionally stable for whistlers; the explicit limit $\Delta t_H\propto h^2$ is 4–8× below the fluid step on the refined meshes. Backward Euler is first order and damps whistlers (dispersion relation error grows like $(\omega\Delta t)^2$ per step). *Spectrum argument (not a proof):* $I+\Delta tMD^T\star_2D=M(M^{-1}+\Delta tL)$, $L=D^T\star_2D$ SPD-semidefinite, is similar to an SPD matrix; the Hall block is approximately skew ($J\times B\cdot J=0$), so the spectrum lies in the right half plane.

**Linear solver.** Restarted GMRES($m=20$) [R22]: Arnoldi with modified Gram–Schmidt, $AV_m=V_{m+1}\bar H_m$, minimise $\|\beta e_1-\bar H_my\|$ by Givens rotations, stop at $\|r\|\le\mathrm{rtol}\,\|E_{\rm expl}+F(\Phi^n)\|$ (true residual recomputed at each restart). Right preconditioning $AP^{-1}u=r_0$, $\delta=P^{-1}V_my$ so the monitored residual is the true one. All vectors are consistent on owners and aliases; inner products run over owned edges + one `Allreduce`.

**Warm start.** The implicit correction $c^n=E^n-E^n_{\rm expl}$ is smooth in time, so the initial guess
$$E^{(0)}_{n+1}=E_{\rm expl}^{\,n+1}+c^{\,n}\ \ (\text{Dirichlet edges: imposed value})$$
has an initial residual $O(\Delta t)\|\dot c\|$ instead of $O(\|c\|)$. *(empirical: sparse mesh 19 → 3.6 it/step, refined mesh 73 → 2.7 it/step, solution identical to the tight one within $4\cdot10^{-5}$.)* This is what makes the implicit step cheaper than the explicit one here.

**Preconditioners.** (a) Jacobi: $P_{ee}=1+\Delta t\,m_e\sum_{f\ni e}\star_2(f)$ (default). (b) Block-Jacobi ILU(0) of the resistive part (optional, `induction_ilu=1`): on each rank the owned-edge matrix $A^{\rm res}_{ee'}=\delta_{ee'}+\Delta t\,m_e\sum_{f\ni e,e'}f2e_{fe}f2e_{fe'}\star_2(f)$ is factored incomplete on its own sparsity pattern, couplings to other ranks and the Hall block dropped (additive Schwarz without overlap [R22]); one extra edge exchange per application refreshes the aliases. Plain ILU(0) [R23] is unstable for large $\Delta t$ because $A^{\rm res}\approx\Delta tML$ is nearly singular on the gradient null space and not an M-matrix; the diagonal is shifted, $\tilde A_{ee}=A_{ee}+\sigma\Delta t\,m_e\,(D^T\star_2D)_{ee}$ ($\sigma=0.3$), the standard remedy for incomplete factorisations of non-M matrices [R24]. *(empirical: ILU(0) 1.7× fewer iterations on the sparse mesh, but not better than Jacobi + warm start on the refined mesh and 10× less accurate at the same rtol; an effective preconditioner for the curl–curl part needs an edge-element multigrid / auxiliary-space method [R25], not implemented.)*

## 10. Verification summary (details in the docs and tests)

| item | check | result |
|---|---|---|
| (1.1), (1.2) | integer products; CT with random EMFs on partitioned/hanging meshes | exact; $|d_2B|/\sum|B_f|\lesssim10^{-15}$–$10^{-12}$ |
| ideal EMF | legacy regression, uniform fields, valence 3–8, np-independence | $4\cdot10^{-15}$; $10^{-14}$; $10^{-15}$ |
| fluid | smooth wave, Sod, Alfvén packet | order 2.0–2.1; $L_1=3.5\cdot10^{-3}$; speed/amplitude/sign correct |
| implicit induction | whistler packet vs. analytic, $\Delta t=77\times\Delta t_{H}$, np=1/4, preconditioners | converges in $h$, $\Delta t$; stable; $10^{-17}$; consistent |
| refinement | topology, B transfer, U conservation | valid; $10^{-14}$–$10^{-11}$; $10^{-14}$ |

## 11. References

Please check the bibliographic details of [R13] and [R14] before citing — they are quoted from memory / from comments in the legacy source; the others are standard and quoted as I know them.

- [R1] A. N. Hirani, *Discrete Exterior Calculus*, PhD thesis, Caltech, 2003.
- [R2] M. Desbrun, A. N. Hirani, M. Leok, J. E. Marsden, *Discrete Exterior Calculus*, arXiv:math/0508341, 2005.
- [R3] R. Hiptmair, *Finite elements in computational electromagnetism*, Acta Numerica 11 (2002) 237–339.
- [R4] C. R. Evans, J. F. Hawley, *Simulation of magnetohydrodynamic flows: a constrained transport method*, Astrophys. J. 332 (1988) 659–677.
- [R5] T. A. Gardiner, J. M. Stone, *An unsplit Godunov method for ideal MHD via constrained transport*, J. Comput. Phys. 205 (2005) 509–539.
- [R6] A. Bossavit, *Computational Electromagnetism: Variational Formulations, Complementarity, Edge Elements*, Academic Press, 1998.
- [R7] C. Burstedde, L. C. Wilcox, O. Ghattas, *p4est: scalable algorithms for parallel adaptive mesh refinement on forests of octrees*, SIAM J. Sci. Comput. 33 (2011) 1103–1133 (octree, 2:1 balance).
- [R8] J. Skilling, *Programming the Hilbert curve*, AIP Conf. Proc. 707 (2004) 381–387.
- [R9] G. Karypis, V. Kumar, *A fast and high quality multilevel scheme for partitioning irregular graphs*, SIAM J. Sci. Comput. 20 (1998) 359–392 (graph model and weighted partitioning; METIS itself is not used in the current build).
- [R10] V. V. Rusanov, *Calculation of interaction of non-steady shock waves with obstacles*, J. Comput. Math. Phys. USSR 1 (1961) 267–279.
- [R11] T. J. Barth, *Recent developments in high order k-exact reconstruction on unstructured meshes*, AIAA Paper 93-0668, 1993.
- [R12] T. J. Barth, D. C. Jespersen, *The design and application of upwind schemes on unstructured meshes*, AIAA Paper 89-0366, 1989.
- [R13] W. Exner et al., 2020, sodium exosphere components (MIV, TD, PSD, SP; Table 2) as cited in the legacy MPCNS_Mercury source comments — **full citation to be confirmed**.
- [R14] Multi-species/multi-fluid MHD for planetary magnetospheres, e.g. X. Jia et al., *Global MHD simulations of Mercury's magnetosphere with coupled planetary interior*, J. Geophys. Res. Space Physics 120 (2015) 4763 — **to be confirmed**; the model used here follows the legacy MPCNS_Mercury formulation.
- [R15] G. Tóth, Y. Ma, T. I. Gombosi, *Hall magnetohydrodynamics on block-adaptive grids*, J. Comput. Phys. 227 (2008) 6967–6984 (Hall term, whistler time-step limit, semi-implicit treatment).
- [R16] D. S. Balsara, D. Spicer, *A staggered mesh algorithm using high order Godunov fluxes to ensure solenoidal magnetic fields in MHD simulations*, J. Comput. Phys. 149 (1999) 270–292.
- [R17] J. D. Huba, *Hall magnetohydrodynamics — a tutorial*, in *Space Plasma Simulation*, Lecture Notes in Physics 615, Springer, 2003, 166–192.
- [R18] S. Gottlieb, C.-W. Shu, E. Tadmor, *Strong stability-preserving high-order time discretization methods*, SIAM Rev. 43 (2001) 89–112.
- [R19] G. Strang, *On the construction and comparison of difference schemes*, SIAM J. Numer. Anal. 5 (1968) 506–517.
- [R20] Linearised/semi-implicit treatment of stiff dispersive terms: see [R15] and the Hall-MHD implicit schemes cited there.
- [R21] D. A. Knoll, D. E. Keyes, *Jacobian-free Newton–Krylov methods: a survey of approaches and applications*, J. Comput. Phys. 193 (2004) 357–397.
- [R22] Y. Saad, M. H. Schultz, *GMRES: a generalized minimal residual algorithm for solving nonsymmetric linear systems*, SIAM J. Sci. Stat. Comput. 7 (1986) 856–869; Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003 (right preconditioning, restarts, Schwarz/ILU).
- [R23] J. A. Meijerink, H. A. van der Vorst, *An iterative solution method for linear systems of which the coefficient matrix is a symmetric M-matrix*, Math. Comp. 31 (1977) 148–162.
- [R24] T. A. Manteuffel, *An incomplete factorization technique for positive definite linear systems*, Math. Comp. 34 (1980) 473–497 (diagonal shift).
- [R25] R. Hiptmair, J. Xu, *Nodal auxiliary space preconditioning in H(curl) and H(div) spaces*, SIAM J. Numer. Anal. 45 (2007) 2483–2509.
