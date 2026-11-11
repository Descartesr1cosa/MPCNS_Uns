#!/usr/bin/env bash
# Build (if needed) and create the run folders of ONE model in a work directory.
#   usage: ./install.sh /path/to/work [Mercury|Mars]     ->  /path/to/work/CASE/01_Meshgen, 02_<Model>Run, ...
# equivalent:  cmake -S . -B build -DUNS_CASE_MODEL=Mars && cmake --build build -j && cmake --install build --prefix /path/to/work
# Existing config.yml / *.txt in the target are kept (re-install only refreshes the programs).
set -e
[ -n "$1" ] || { echo "usage: $0 <work-dir> [Mercury|Mars]"; exit 2; }
ROOT=$(cd "$(dirname "$0")" && pwd); WORK=$(mkdir -p "$1" && cd "$1" && pwd)
BUILD=${BUILD:-$ROOT/build}
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DUNS_BUILD_TESTS=${TESTS:-OFF} -DUNS_CASE_MODEL="${2:-Mercury}"
cmake --build "$BUILD" -j "${JOBS:-16}"
cmake --install "$BUILD" --prefix "$WORK"
