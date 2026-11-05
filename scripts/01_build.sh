#!/usr/bin/env bash
# Configure and build (C++20, MPI).  Usage: scripts/01_build.sh
set -e
source "$(dirname "$0")/env.sh"
cmake -S "$UNS_ROOT" -B "$UNS_BUILD" -DCMAKE_BUILD_TYPE=Release
cmake --build "$UNS_BUILD" -j "${JOBS:-16}"
echo "built: $UNS_BUILD/Preprocess/mpcns_pre  $UNS_BUILD/MercuryUns/MercuryUns"
