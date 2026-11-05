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
} // namespace pre
