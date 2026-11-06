#pragma once
#include "uns/parallel/LocalMesh.h"

namespace mhd
{
struct DivBReport
{
    double max_abs = 0;   // max_c |(d2 B)_c|  (flux imbalance of a cell)
    double max_rel = 0;   // relative to sum_f |B_f| of the cell
    long worst_cell_gid = -1;
};
// Over the owned cells of this rank (global max is the caller's MPI_Allreduce).
// `scale_extra`: optional second flux array (background field) added to the normalisation sum|B_f| + |B_extra,f|, so that the
// round-off level is measured against the total field the update works with.
DivBReport div_b(const uns::CellComplex &cx, const double *B, const double *scale_extra = nullptr);
} // namespace mhd
