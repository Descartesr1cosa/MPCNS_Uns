#pragma once
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
struct DivBReport
{
    double max_abs = 0;   // max_c |(d2 B)_c|  (flux imbalance of a cell)
    double max_rel = 0;   // relative to sum_f |B_f| of the cell
    long worst_cell_gid = -1;
};
// Over the owned cells of this rank (global max is the caller's MPI_Allreduce).
DivBReport div_b(const uns::CellComplex &cx, const double *B);
} // namespace mercury
