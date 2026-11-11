#!/usr/bin/env bash
# End-to-end smoke test: generate a partitioned cubed-sphere shell, then run the MercuryUns
# CT self-test on it with several process counts.  The printed checksum(B) must be identical.
#   usage: tools/run_selftest.sh [build-dir] [work-dir]
set -e
B=${1:-build}; W=${2:-${TMPDIR:-/tmp}/uns_selftest}
mkdir -p "$W"
for n in 1 2 5 8; do
  "$B/Preprocess/mpcns_pre" tools/demo/demo_mesh.txt case_template/common/partition.txt output_dir="$W/g$n" proc_num=$n write_vtk=0 partition_method=rcb > /dev/null
  mpirun --allow-run-as-root --oversubscribe -np $n "$B/MercuryUns/MercuryUns" - geometry_dir="$W/g$n" | tail -1
done
