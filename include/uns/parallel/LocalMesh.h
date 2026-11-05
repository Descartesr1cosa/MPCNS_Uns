#pragma once
// Per-rank sub-complex plus its communication plan.
//
// Local numbering of every dimension is [owned | ghost]; owned entities are those
// whose owner is this rank (owner = lowest rank among the ranks of the cells
// incident to the entity), ghosts are sorted by (owner, global id).  The local
// cells are the owned cells plus one vertex-neighbour layer of ghost cells, so
// every edge that touches an owned cell has *all* its incident cells/faces locally.
//
// Plans (all index lists are local ids, ordered by ascending global id so that
// sender and receiver agree without any negotiation):
//   send[d] / recv[d] : owner -> alias copy ("forward"), per dimension
//   reduce_send/recv  : edges only; alias -> owner partial sums (edge assembly)
#include <array>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns
{
struct CommLists
{
    std::vector<int> peers;
    std::vector<std::int64_t> ptr{0};
    std::vector<Id> idx;
    std::size_t total() const { return idx.size(); }
    void add_peer(int peer, const std::vector<Id> &ids)
    {
        peers.push_back(peer);
        idx.insert(idx.end(), ids.begin(), ids.end());
        ptr.push_back(static_cast<std::int64_t>(idx.size()));
    }
};

struct CommPlan
{
    std::array<CommLists, kNumDims> send, recv;
    CommLists reduce_send, reduce_recv;
};

struct LocalMesh
{
    Mesh mesh;
    int rank = 0, nranks = 1;
    CommPlan plan;
    Id n_owned(Dim d) const { return mesh.cx.ent[di(d)].n_owned; }
    Id n_total(Dim d) const { return mesh.cx.count(d); }
};
} // namespace uns
