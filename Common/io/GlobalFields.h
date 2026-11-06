#pragma once
// Gather owned-entity data to rank 0 in the global numbering (coarse solution / restart files).
#include <mpi.h>

#include <vector>

#include "uns/complex/CellComplex.h"

namespace mhd
{
// values: ncomp doubles per local entity (owned first).  Rank 0 returns n_global*ncomp values, the others an empty vector.
std::vector<double> gather_entities(MPI_Comm comm, const uns::EntityInfo &ent, const double *values, int ncomp);
} // namespace mhd
