# CASE — how to work with it

Created by `./install.sh <work-dir> [Mercury|Mars]` (or `cmake --install build --prefix <work-dir>`): one model per CASE.
Every step is a numbered folder: **enter it, check `config.yml` and the `*.txt`, run the one program in it.**

```
Mercury:  01_Meshgen ./mpcns_pre  ->  02_MercuryRun ./MercuryUns  ->  03_MercuryRefine ./mpcns_pre  ->  04_MercuryRunFine ./MercuryUns
Mars:     01_Meshgen ./mpcns_pre  ->  02_MarsRun ./MarsUns
```

## Where does a parameter live?  (every key exists in exactly one file)

`config.yml` says **where** and **how many ranks**; the txt files say **what and how**, one file per kind of setting:

| file | kind of setting | examples |
|---|---|---|
| `config.yml` | paths, task, parallel launch, which txt files to read | `geometry_dir`, `output_dir`, `init_fields`, `nproc`, `mpi_launcher`, `Crustal_Dipoles_File`, `settings:` |
| `mesh.txt` (01) | mesh generation method and sizes | `cs_ncube`, `cs_nr`, `cs_r_stretch` (growth rate), `cs_ell_lo/hi` |
| `partition.txt` (01) | partition method, cost model (03 re-uses 01's) | `partition_method`, `cost_*` |
| `refine.txt` (03) | refinement regions and indicators | `refine_max_level`, `refine_region_N`, `refine_indicator_N` |
| `physics.txt` (02) | constants, solar wind, planet, model switches, boundary conditions | `sw_*`, `dipole_*`, `is_Hall_Efield`, `bc_*` |
| `numerics.txt` (02) | time step, order, limiter, induction solver, steady-state options | `CFL`, `second_order`, `induction_rtol` |
| `run.txt` (02, 04) | end time, logging and output | `max_Time`, `log_every`, `output_step`, `restart_step` |

Step 04 reads the physics and numerics of step 02 (single copy) and has its own `run.txt`; step 03 reads 01's `partition.txt`.
`config.yml` format: flat `key: value`, `#` comments; `settings: [a.txt, b.txt]` lists the txt files read (later wins, then the yml keys).
txt format: `type name value  // comment`.

## Notes
* `./<program>` reads `./config.yml` (or `./<program> other.yml`); relative paths are relative to the folder.  Any key can be overridden on the command line: `./MercuryUns max_Time=5`.
* The solver starts MPI itself: `nproc: auto` = number of partitions found in `geometry_dir`.  Inside a batch job (`srun`/`mpirun` already active) it does not relaunch itself.
* Runtime dependency: only an MPI library (Linux: the one used at build time, found through the embedded RPATH or `mpirun` in `PATH`; Windows: MS-MPI).  All other code is linked statically.
* Re-running the install keeps the config.yml / txt you edited and only refreshes the programs.
* Restart: set `init_fields:` to a `fields_final.unsb` of an earlier run on the same mesh.  Results: `DATA/solution.pvd` (ParaView).

## Windows
* **Recommended: WSL2** (Ubuntu): `sudo apt install build-essential cmake openmpi-bin libopenmpi-dev`, then exactly the Linux steps (`./install.sh ~/work Mercury`).  Keep the CASE inside the WSL file system (fast I/O); Windows reaches it as `\\wsl$\Ubuntu\home\<you>\work\CASE`, ParaView for Windows opens `DATA/solution.pvd` from there.
* **Native Windows** (MSVC or MSYS2/MinGW + MS-MPI SDK): `install.cmd C:\work Mercury` from a Developer Command Prompt (same as the cmake commands in `install.sh`).  The code uses only portable C++20 (`std::filesystem`), the self-launch uses `mpiexec` (MS-MPI) and runs `MercuryUns.exe` as the launcher's child; the C/C++ runtime is linked statically, so only the MS-MPI runtime (`msmpi.dll`, MS-MPI redistributable) is needed on the target PC.  In `config.yml` use `/` in paths (`../01_Meshgen/geometry`).  Start with `.\MercuryUns.exe` (PowerShell) or `MercuryUns.exe`.  **Not tested on Windows** (no Windows toolchain was available when this was written): please report build errors.
