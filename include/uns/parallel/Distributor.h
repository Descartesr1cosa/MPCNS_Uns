#pragma once
#include <vector>

#include "uns/parallel/LocalMesh.h"

namespace uns
{
// Splits a global Mesh according to a cell partition into per-rank LocalMeshes.
class Distributor
{
public:
    Distributor(const Mesh &global, const std::vector<int> &cell_rank, int nranks);
    LocalMesh build(int rank) const;
    // owner rank of every global entity of dimension d
    const std::vector<int> &owner(Dim d) const { return owner_[di(d)]; }

private:
    const Mesh &g_;
    int nranks_;
    std::array<std::vector<int>, kNumDims> owner_;
    // per rank, per dim: local entity list = [owned (gid order) | ghost (owner,gid order)]
    std::vector<std::array<std::vector<GId>, kNumDims>> local_;
    std::vector<std::array<Id, kNumDims>> n_owned_;
    // edges touched by owned cells of a rank but owned elsewhere, sorted (owner, gid)
    std::vector<std::vector<GId>> touch_edges_;
};
} // namespace uns
