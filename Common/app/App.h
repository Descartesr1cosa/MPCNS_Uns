#pragma once
// Common driver of the solver executables (MercuryUns, MarsUns):
//   <exe> setup-file... [name=value ...]     one MPI rank per part_XXXXX.unsb written by mpcns_pre
//   task=selftest   structural test of the induction path (B -= dt d1 E keeps d2 B at round-off)
//   task=run        calls the model's run function
#include <mpi.h>

#include <functional>

#include "uns/core/Params.h"
#include "uns/parallel/LocalMesh.h"

namespace mhd
{
using RunFn = std::function<void(MPI_Comm, const uns::LocalMesh &, const uns::Params &)>;
int run_app(int argc, char **argv, const char *name, const RunFn &run);
} // namespace mhd
