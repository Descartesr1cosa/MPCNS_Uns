#!/usr/bin/env bash
# Run on the refined mesh (first order in space, semi-implicit induction with Hall + resistivity + ambipolar field).
#   usage: scripts/04_fine_run.sh [NP=56] [T_END=50]     (NP must equal the NP of 03_refine.sh)
# Restart: set RESTART=1 to continue from $CASE_DIR/fine/DATA/fields_final.unsb (rewritten at the end of every run and every RESTART_EVERY steps).
set -e
source "$(dirname "$0")/env.sh"
NP=${1:-56}; TEND=${2:-50}
F=$CASE_DIR/fine; SOL=$UNS_BUILD/MercuryUns/MercuryUns; IN=$UNS_ROOT
INIT=$F/geometry/fine_fields.unsb
if [ "${RESTART:-0}" = "1" ] && [ -f $F/DATA/fields_final.unsb ]; then cp $F/DATA/fields_final.unsb $F/restart.unsb; INIT=$F/restart.unsb; fi
cd $F
$MPIRUN $NP $SOL $IN/MercuryUns/INPUT/setup.txt $IN/MercuryUns/INPUT/mercury.txt \
    geometry_dir=$F/geometry init_fields=$INIT output_dir=$F/DATA max_Time=$TEND second_order=0 \
    log_every=${LOG_EVERY:-1000} output_step=${OUT_EVERY:-20000} restart_step=${RESTART_EVERY:-20000} | tee -a $F/run.log
