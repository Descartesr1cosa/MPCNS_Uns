# Common settings for the Mercury workflow scripts.  Edit the paths below for your cluster.
# MPI / compilers: load your modules here, e.g.  module load gcc/13 openmpi/4.1.6
export UNS_ROOT=${UNS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
export UNS_BUILD=${UNS_BUILD:-$UNS_ROOT/build}
# case directory with the working data (coarse mesh, coarse solution, refined mesh, outputs)
export CASE_DIR=${CASE_DIR:-$PWD/mercury_case}
# MPI launcher.  Examples:  "mpirun -np" (default) | "mpirun --map-by core -np" | "srun -n"
export MPIRUN=${MPIRUN:-"mpirun -np"}
