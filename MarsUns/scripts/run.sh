#!/usr/bin/env bash
# MarsUns: coarse mesh without solid shell + ten-ion run.   usage: MarsUns/scripts/run.sh [NP=4] [T_END=0.1] [extra name=value ...]
# e.g.  MarsUns/scripts/run.sh 8 1.0 is_Crustal_Bfield=1 Crustal_Dipoles_File=$PWD/MarsUns/INPUT/Dipoles_4840
set -e
source "$(dirname "$0")/../../scripts/env.sh"
NP=${1:-4}; TEND=${2:-0.1}; shift 2 2>/dev/null || true
C=${MARS_DIR:-$UNS_ROOT/Run/Mars}; mkdir -p "$C"; cd "$C"
$UNS_BUILD/Preprocess/mpcns_pre $UNS_ROOT/Preprocess/INPUT/setup.txt $UNS_ROOT/MarsUns/INPUT/mars_mesh.txt output_dir=$C/geometry proc_num=$NP write_vtk=0
$MPIRUN $NP $UNS_BUILD/MarsUns/MarsUns $UNS_ROOT/MarsUns/INPUT/mars.txt \
    geometry_dir=$C/geometry output_dir=$C/DATA max_Time=$TEND "$@" | tee $C/run.log
