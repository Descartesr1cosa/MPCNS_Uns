#!/usr/bin/env bash
# Refine the coarse mesh from the coarse flow field and partition it for NP ranks with the measured cost model.
#   usage: scripts/03_refine.sh [NP=56]
# The refinement regions / indicators below are the ones used for the ~0.57M-cell mesh (first fluid layer 0.0063 R_M); edit freely.
set -e
source "$(dirname "$0")/env.sh"
NP=${1:-56}
C=$CASE_DIR/coarse; F=$CASE_DIR/fine; mkdir -p "$F"
PRE=$UNS_BUILD/Preprocess/mpcns_pre; IN=$UNS_ROOT
$PRE $IN/Preprocess/INPUT/setup.txt task=refine coarse_dir=$C/geometry coarse_fields=$C/DATA/fields_final.unsb output_dir=$F/geometry proc_num=$NP write_vtk=1 \
    cost_profile=$C/DATA/load_profile.txt \
    refine_max_level=3 refine_solid_max_level=0 refine_spherical_snap=1 \
    "refine_region_1=shell 1.0 1.05 3" "refine_region_2=shell 1.05 1.15 2" "refine_region_3=shell 1.15 1.4 1" \
    "refine_region_4=box -5 -1.6 -1.6 -1.2 1.6 1.6 1" \
    "refine_indicator_1=current 0.5 2" "refine_indicator_2=density 0.3 1" | tee $F/refine.log
echo "refined mesh: $F/geometry (partition.pvtu shows the ranks)"
