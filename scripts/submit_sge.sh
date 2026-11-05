#!/bin/bash
#$ -N mercury_uns
#$ -cwd
#$ -j y
#$ -pe mpi 56
#$ -q normal.q
# SGE job script (like the legacy runs: NSLOTS ranks).  Submit with:   qsub scripts/submit_sge.sh
# Adjust the queue (-q), the parallel environment (-pe) and the module lines to your cluster.
# module load gcc/13.2.0 openmpi/4.1.6
export CASE_DIR=${CASE_DIR:-$PWD/mercury_case}
source "$(dirname "$0")/env.sh"
NP=${NSLOTS:-56}
# first time:  scripts/01_build.sh ; scripts/02_coarse_run.sh 28 25 ; scripts/03_refine.sh $NP   (or submit those as separate smaller jobs)
RESTART=${RESTART:-0} scripts/04_fine_run.sh $NP ${T_END:-50}
