# MPCNS_Uns
Fully unstructured, DEC-preserving Hall-MHD on a signed cell complex (cell -> face -> edge -> node), for the magnetospheres of
Mercury (two ions, H+ / Na+) and Mars (ten ions).  See [PLAN.md](PLAN.md) for the plan, design decisions and status.

## Quick start: install a ready-to-run CASE
```
./install.sh /path/to/work Mercury     # or Mars;  = cmake build + cmake --install build --prefix /path/to/work -DUNS_CASE_MODEL=...
cd /path/to/work/CASE/01_Meshgen && ./mpcns_pre          # then 02_MercuryRun: ./MercuryUns, 03_MercuryRefine, 04_MercuryRunFine
```
Each step folder holds the program, a `config.yml` (paths, ranks) and the parameter files by kind (`mesh/partition/refine/physics/numerics/run.txt`); the solver starts MPI itself.  Which parameter lives where: [case_template/README.md](case_template/README.md) (templates = the only copy of the input files).

## Windows
* **Recommended: WSL2** (Ubuntu): `sudo apt install build-essential cmake openmpi-bin libopenmpi-dev`, then exactly the Linux steps (`./install.sh ~/work Mercury`).  Keep the CASE inside the WSL file system (fast I/O); Windows reaches it as `\\wsl$\Ubuntu\home\<you>\work\CASE`, ParaView for Windows opens `DATA/solution.pvd` from there.
* **Native Windows** (MSVC or MSYS2/MinGW + MS-MPI SDK): `install.cmd C:\work Mercury` from a Developer Command Prompt (same as the cmake commands in `install.sh`).  The code uses only portable C++20 (`std::filesystem`), the self-launch uses `mpiexec` (MS-MPI) and runs `MercuryUns.exe` as the launcher's child; the C/C++ runtime is linked statically, so only the MS-MPI runtime (`msmpi.dll`, MS-MPI redistributable) is needed on the target PC.  In `config.yml` use `/` in paths (`../01_Meshgen/geometry`).  Start with `.\MercuryUns.exe` (PowerShell) or `MercuryUns.exe`.  **Not tested on Windows** (no Windows toolchain was available when this was written): please report build errors.

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
                two-ion ambipolar term, interior resistivity
MarsUns/        Mars only (namespace mars): driver/MarsSolver, physics/ (ten ions, photochemistry, atmosphere, crustal dipoles),
                docs/, tools/, tests/
tests/          library + Mercury + shared-block tests;   MarsUns/tests/  Mars tests
docs/           Mercury-side and general documentation;   MarsUns/docs/Mars.md  Mars branch
case_template/  input templates installed by install.sh (config.yml + txt per step)
```
Executables: `build/Preprocess/mpcns_pre`, `build/MercuryUns/MercuryUns`, `build/MarsUns/MarsUns`.

Docs: [PLAN.md](PLAN.md) (plan + milestone status), docs/MercuryRun.md (Mercury state, open problems, read first when resuming), docs/Methods.md (math), docs/Multigrid.md, docs/Solver.md, MarsUns/docs/Mars.md.
Tests: `cd build && ctest` (23 tests).

Output: ParaView reads the VTU/PVD files directly; for Tecplot convert with `tools/vtu2tec/vtu2tec.py` driven by `tools/vtu2tec/config.yml` (binary .plt, see docs/FileFormats.md).
