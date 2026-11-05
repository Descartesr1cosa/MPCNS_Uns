#!/usr/bin/env bash
# Sparse mesh + first-order coarse run.  Writes the coarse mesh (global.unsb), the coarse solution (fields_final.unsb) and load_profile.json/.txt.
#   usage: scripts/02_coarse_run.sh [NP=28] [T_END=25]
set -e
source "$(dirname "$0")/env.sh"
NP=${1:-28}; TEND=${2:-25}
C=$CASE_DIR/coarse; mkdir -p "$C"; cd "$C"
PRE=$UNS_BUILD/Preprocess/mpcns_pre; SOL=$UNS_BUILD/MercuryUns/MercuryUns; IN=$UNS_ROOT
$PRE $IN/Preprocess/INPUT/setup.txt $IN/Preprocess/INPUT/mercury_mesh.txt output_dir=$C/geometry proc_num=$NP cs_r_stretch=1.07 write_vtk=0
$MPIRUN $NP $SOL $IN/MercuryUns/INPUT/setup.txt $IN/MercuryUns/INPUT/mercury.txt \
    geometry_dir=$C/geometry output_dir=$C/DATA max_Time=$TEND second_order=0 log_every=2000 output_step=0 | tee $C/run.log
echo "coarse solution: $C/DATA/fields_final.unsb ; load profile: $C/DATA/load_profile.txt"
