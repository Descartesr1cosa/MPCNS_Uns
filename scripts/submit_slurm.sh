#!/bin/bash
#SBATCH -J mercury_uns
#SBATCH -n 56
#SBATCH -t 48:00:00
#SBATCH -o mercury_uns.%j.log
# SLURM job script.  Submit with:   sbatch scripts/submit_slurm.sh
# module load gcc/13.2.0 openmpi/4.1.6
export CASE_DIR=${CASE_DIR:-$PWD/mercury_case}
export MPIRUN="srun -n"
source "$(dirname "$0")/env.sh"
RESTART=${RESTART:-0} scripts/04_fine_run.sh ${SLURM_NTASKS:-56} ${T_END:-50}
