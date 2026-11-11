#!/usr/bin/env bash
# Coarse shell -> synthetic coarse solution -> refinement -> partition -> CT self-test on the refined mesh.
#   usage: tools/run_refine_demo.sh [build-dir] [work-dir]
set -e
B=${1:-build}; W=${2:-${TMPDIR:-/tmp}/uns_refine_demo}
mkdir -p "$W"
PRE="$B/Preprocess/mpcns_pre"
$PRE tools/demo/demo_mesh.txt case_template/common/partition.txt output_dir="$W/coarse/geometry" proc_num=4 cs_ncube=8 cs_nr=10 cs_nsolid=2 cs_r_in=0.8 cs_r_out=14 | grep -E "topology|geometry"
$PRE - task=demo_fields coarse_dir="$W/coarse/geometry"
$PRE tools/demo/demo_mesh.txt case_template/common/partition.txt tools/demo/demo_refine.txt task=refine coarse_dir="$W/coarse/geometry" \
     coarse_fields="$W/coarse/geometry/coarse_fields.unsb" output_dir="$W/fine/geometry" proc_num=8 | grep -E "refine\]|level |transfer|imbalance"
mpirun --allow-run-as-root --oversubscribe -np 8 "$B/MercuryUns/MercuryUns" - geometry_dir="$W/fine/geometry" | tail -1
echo "ParaView: open $W/fine/geometry/partition.pvtu (colour by level / rank / rho)"
