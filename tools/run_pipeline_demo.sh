#!/usr/bin/env bash
# Whole workflow on a small Mercury-like case (all rank counts <= 84):
#   coarse complex -> coarse MHD run -> refinement from the coarse solution -> run on the refined non-conforming mesh
#   usage: tools/run_pipeline_demo.sh [build-dir] [work-dir]
set -e
B=$(cd "${1:-build}" && pwd); W=${2:-${TMPDIR:-/tmp}/uns_pipeline}; mkdir -p "$W"; W=$(cd "$W" && pwd)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PRE="$B/Preprocess/mpcns_pre"; SOLVER="$B/MercuryUns/MercuryUns"; SETUP="$ROOT/case_template/Mercury/02_MercuryRun/physics.txt $ROOT/case_template/Mercury/02_MercuryRun/numerics.txt $ROOT/case_template/Mercury/02_MercuryRun/run.txt task=run second_order=1"
SHELL_ARGS="cs_ncube=6 cs_nr=8 cs_nsolid=2 cs_r_in=0.7 cs_r_out=14 cs_r_stretch=1.15"

echo "== 1. coarse mesh (4 parts)"
$PRE "$ROOT/tools/demo/demo_mesh.txt" "$ROOT/case_template/common/partition.txt" output_dir="$W/coarse/geometry" proc_num=4 $SHELL_ARGS write_vtk=0 | grep -E "topology|imbalance"
echo "== 2. coarse run"
( cd "$W/coarse" && mpirun --allow-run-as-root --oversubscribe -np 4 $SOLVER $SETUP geometry_dir="$W/coarse/geometry" output_dir="$W/coarse/DATA" max_Time=${T_COARSE:-6} log_every=40 | tail -3 )
echo "== 3. refine (flow indicators + tail region), 8 parts"
$PRE "$ROOT/tools/demo/demo_mesh.txt" "$ROOT/case_template/common/partition.txt" "$ROOT/tools/demo/demo_refine.txt" task=refine coarse_dir="$W/coarse/geometry" \
     coarse_fields="$W/coarse/DATA/fields_final.unsb" output_dir="$W/fine/geometry" proc_num=8 \
     refine_indicator_1="current 0.5 2" refine_indicator_2="density 0.3 1" | grep -E "refine\]|level |transfer|imbalance"
echo "== 4. run on the refined mesh (restart from the transferred coarse solution)"
( cd "$W/fine" && mpirun --allow-run-as-root --oversubscribe -np 8 $SOLVER $SETUP geometry_dir="$W/fine/geometry" output_dir="$W/fine/DATA" \
     init_fields="$W/fine/geometry/fine_fields.unsb" max_Time=${T_FINE:-8} log_every=40 | tail -4 )
echo "ParaView: $W/fine/DATA/solution.pvd  and  $W/fine/geometry/partition.pvtu"
