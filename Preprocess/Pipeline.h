#pragma once
#include "uns/core/Params.h"

namespace pre
{
// task = generate : build/read a mesh, partition it, write per-rank sub-complexes + report
// task = inspect  : validate a complex written by `generate` (topology + geometry)
// task = export   : write ParaView files (.vtu/.pvtu/.pvd) for an existing case directory
int run_generate(const uns::Params &p);
int run_inspect(const uns::Params &p);
int run_export(const uns::Params &p);
// task = refine      : coarse complex (+ coarse solution) -> refined complex, transferred B/U, new partition
// task = demo_fields : synthetic coarse solution for exercising `refine` (until the solver writes real ones)
int run_refine(const uns::Params &p);
int run_demo_fields(const uns::Params &p);
} // namespace pre
