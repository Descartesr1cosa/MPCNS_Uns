#pragma once
// Reading/writing meshes, sub-complexes and entity-numbered fields as UNSBIN files.
#include <string>
#include <vector>

#include "uns/parallel/LocalMesh.h"

namespace uns::io
{
// A Mesh is stored as a LocalMesh with nranks = 1 and empty plans.
void write_local_mesh(const std::string &path, const LocalMesh &lm);
LocalMesh read_local_mesh(const std::string &path);

// Field data attached to the entities of a complex (restart / coarse solution).
// `values` has count(dim) * ncomp doubles in the *local* numbering of the file it belongs to.
struct FieldRecord
{
    std::string name;
    Dim dim = Dim::Cell;
    int ncomp = 1;
    std::vector<double> values;
};
// `topology_hash` ties the field file to the complex it was computed on.
void write_fields(const std::string &path, std::uint64_t topology_hash, const std::vector<FieldRecord> &fields);
std::vector<FieldRecord> read_fields(const std::string &path, std::uint64_t *topology_hash = nullptr);
} // namespace uns::io
