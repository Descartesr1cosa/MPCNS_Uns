#pragma once
// ParaView/VisIt-readable XML unstructured grid output (binary appended raw data).
// Every cell is written as VTK_POLYHEDRON (type 42), so hanging-node polyhedra,
// hexahedra and future cell types need no special casing.
#include <string>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns::io
{
struct VtkField
{
    std::string name;
    int ncomp = 1;
    std::vector<double> values; // one tuple per cell (or per node)
};

// Writes the first `n_cells_out` cells (the owned ones) of `m`.
void write_vtu(const std::string &path, const Mesh &m, Id n_cells_out, const std::vector<VtkField> &cell_data,
               const std::vector<VtkField> &point_data = {});
// Parallel master file referencing per-rank pieces (paths relative to the master).
void write_pvtu(const std::string &path, const std::vector<std::string> &pieces, const std::vector<VtkField> &cell_data,
                const std::vector<VtkField> &point_data = {});
// Time-series collection.
void write_pvd(const std::string &path, const std::vector<std::pair<double, std::string>> &entries);
} // namespace uns::io
