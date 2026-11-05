#!/usr/bin/env bash
# Steady-state (pseudo-time) run with local time stepping and the FAS multigrid correction (docs/Multigrid.md).
#   usage: scripts/05_steady_run.sh [STAGE=coarse|fine] [NP] [MAX_STEPS]
#   coarse : mesh of 02_coarse_run.sh (partitioned for NP, default 12), starts from the analytic initial state   -> $CASE_DIR/coarse_steady/DATA/fields_final.unsb
#   fine   : mesh of 03_refine.sh (NP must be its NP), starts from its fine_fields.unsb (or RESTART=1 from fine_steady/DATA)
# MODE=fast (default): B_RATIO=32 + 4-level multigrid: ~20x fewer iterations than the time-accurate run, but the converged state is 15-70 % off (docs/Multigrid.md 6.5):
#                      use it to get a first approximation (e.g. the coarse solution for the refinement indicators).
# MODE=accurate:       B_RATIO=4, no multigrid (it does not help once B limits the convergence): ~8x fewer iterations, a few percent error in rho, B.
# Other knobs (environment): STEADY_TOL (stop when res_H, res_B < tol, default 1e-7), MG_LEVELS, MG_EVERY (4), MG_START (200), MG_RELAX (1; 0.5 removes a residual
#                      floor of Na+ at small B_RATIO), B_RATIO (cap of the induction step in units of the smallest cell step; > 32 destabilises), LOG_EVERY (100)
# For the exact steady state use the time-accurate run (04_fine_run.sh, ~24000 steps on the coarse mesh).
set -e
source "$(dirname "$0")/env.sh"
STAGE=${1:-coarse}; NP=${2:-12}; MAXSTEP=${3:-5000}
PRE=$UNS_BUILD/Preprocess/mpcns_pre; SOL=$UNS_BUILD/MercuryUns/MercuryUns; IN=$UNS_ROOT
MODE=${MODE:-fast}
if [ "$MODE" = accurate ]; then DEF_R=4; DEF_L=1; else DEF_R=32; DEF_L=4; fi
OPTS="steady=1 second_order=0 lts_b_ratio=${B_RATIO:-$DEF_R} mg_levels=${MG_LEVELS:-$DEF_L} mg_cycle_every=${MG_EVERY:-4} mg_start=${MG_START:-200} mg_relax=${MG_RELAX:-1} steady_tol=${STEADY_TOL:-1e-7} max_Nstep=$MAXSTEP log_every=${LOG_EVERY:-100} write_vtk=0"
if [ "$STAGE" = coarse ]; then
  D=$CASE_DIR/coarse_steady; mkdir -p "$D"; cd "$D"
  $PRE $IN/Preprocess/INPUT/setup.txt $IN/Preprocess/INPUT/mercury_mesh.txt output_dir=$D/geometry proc_num=$NP cs_r_stretch=1.07 write_vtk=0
  $MPIRUN $NP $SOL $IN/MercuryUns/INPUT/setup.txt $IN/MercuryUns/INPUT/mercury.txt geometry_dir=$D/geometry output_dir=$D/DATA $OPTS | tee $D/run.log
else
  D=$CASE_DIR/fine; cd "$D"; mkdir -p steady
  INIT=$D/geometry/fine_fields.unsb
  if [ "${RESTART:-0}" = "1" ] && [ -f $D/steady/DATA/fields_final.unsb ]; then cp $D/steady/DATA/fields_final.unsb $D/steady/restart.unsb; INIT=$D/steady/restart.unsb; fi
  $MPIRUN $NP $SOL $IN/MercuryUns/INPUT/setup.txt $IN/MercuryUns/INPUT/mercury.txt geometry_dir=$D/geometry init_fields=$INIT output_dir=$D/steady/DATA $OPTS | tee -a $D/steady/run.log
fi
echo "convergence table:  python3 $IN/tools/conv_report.py <run.log>"
