# Mercury workflow scripts (run them yourself)

All scripts read `scripts/env.sh` (paths, MPI launcher `MPIRUN`, default `mpirun -np`; use `MPIRUN="srun -n"` under SLURM). The working data go to `$CASE_DIR` (default `./mercury_case`).

| step | command | what it does | time (this 256-core node) |
|---|---|---|---|
| 1 | `scripts/01_build.sh` | CMake Release build (needs C++20, MPI) | ~2 min |
| 2 | `scripts/02_coarse_run.sh 28 25` | sparse mesh (24192 cells), first order to t = 25 on 28 ranks; writes `coarse/DATA/fields_final.unsb` and `load_profile.json/.txt` | ~5 min |
| 3 | `scripts/03_refine.sh 56` | refinement from the coarse flow field (574665 cells, 45924 polyhedra, first fluid layer 0.0063 R_M), partition for 56 ranks with the measured cost model | ~1 min |
| 4 | `scripts/04_fine_run.sh 56 50` | run on the refined mesh to t = 50; `RESTART=1` continues from `fine/DATA/fields_final.unsb`; knobs: `OUT_EVERY` (ParaView dumps), `RESTART_EVERY`, `LOG_EVERY` | 87 ms/step, dt ~ 8e-5 -> ~7 h per 25 time units |
| 5 | `scripts/05_steady_run.sh coarse 12` / `fine 56 20000` | steady-state (pseudo-time) run with local steps + FAS multigrid, stops at `STEADY_TOL`; see `docs/Multigrid.md` | coarse mesh: ~45 s on 3 ranks to res 1e-5 |
Job scripts: `scripts/submit_sge.sh` (`qsub`, like the legacy `NSLOTS` runs) and `scripts/submit_slurm.sh` (`sbatch`); adjust queue / modules. The number of ranks of step 4 must equal that of step 3.
Outputs: `fine/DATA/solution.pvd` (ParaView), `fields_final.unsb` (restart), `load_profile.json` (timings per rank), `run.log`.
