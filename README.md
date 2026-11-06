# MPCNS_Uns
Fully unstructured, DEC-preserving Hall-MHD on a signed cell complex (cell -> face -> edge -> node), for the magnetospheres of
Mercury (two ions, H+ / Na+) and Mars (ten ions).  See [PLAN.md](PLAN.md) for the plan, design decisions and status.

## Layout (dependencies point downwards only)
```
src/uns/        mesh library: core, complex, geometry, generate, partition, parallel, refine, io, comm   (namespace uns)
Preprocess/     mpcns_pre: mesh generation / partition / refinement tool
Common/         model-independent MHD blocks shared by both solvers                                     (namespace mhd, library mhd_common)
  config/         Config: run configuration read from setup files
  physics/        Euler helpers + Rusanov flux, nondimensional scales, chemistry-rate struct
  state/          FieldSet
  boundary/       face / edge classification (wall, inflow, outflow)
  fluid/          limited-gradient reconstruction, finite-volume flux right-hand side
  induction/      constrained transport (B -= dt d1 E), edge assembly, ideal UCT EMF, cell current,
                  semi-implicit Hall / resistive solve (GMRES), Hall hook
  diagnostics/    div B            io/  global-field gather (restart)            app/  common main(): setup, mesh read, selftest
MercuryUns/     Mercury only (namespace mercury): driver/Solver, multigrid, initial state, Na+ exosphere/sources, gyration,
                two-ion ambipolar term, interior resistivity;  INPUT/ (setup, mercury.txt, mercury_mesh.txt)
MarsUns/        Mars only (namespace mars): driver/MarsSolver, physics/ (ten ions, photochemistry, atmosphere, crustal dipoles),
                INPUT/, docs/, scripts/, tools/, tests/
tests/          library + Mercury + shared-block tests;   MarsUns/tests/  Mars tests
docs/           Mercury-side and general documentation;   MarsUns/docs/Mars.md  Mars branch
scripts/        Mercury workflow (01-05);                 MarsUns/scripts/run.sh
```
Executables: `build/Preprocess/mpcns_pre`, `build/MercuryUns/MercuryUns`, `build/MarsUns/MarsUns`.

Output: ParaView reads the VTU/PVD files directly; for Tecplot convert with `tools/vtu2tec/vtu2tec.py` driven by `tools/vtu2tec/config.yml` (binary .plt, see docs/FileFormats.md).
