// partition_cells_grouped: groups are never split, the loads of the ranks are balanced on the quotient graph
#include <cstdio>

#include "TestUtil.h"
#include "uns/generate/HexMesh.h"
#include "uns/partition/Partitioner.h"

using namespace uns;

int main()
{
    BoxSpec bs;
    bs.nx = 16; bs.ny = 16; bs.nz = 8; bs.hi = {16, 16, 8};
    const Mesh m = build_mesh(make_box(bs));
    const Geometry g = compute_geometry(m);
    const CellGraph cg = build_cell_graph(m, g, CostModel());
    // 4x4x4 blocks of the box (cell order of make_box is not assumed: blocks from the cell centres)
    std::vector<Id> grp(m.cx.n_cells());
    for (Id c = 0; c < m.cx.n_cells(); ++c)
    {
        const Vec3 x = g.cell_center[c];
        grp[c] = (static_cast<Id>(x.x / 4) * 4 + static_cast<Id>(x.y / 4)) * 2 + static_cast<Id>(x.z / 4);
    }
    PartOptions po;
    for (int np : {3, 5, 8})
    {
        const std::vector<int> part = partition_cells_grouped(cg, grp, np, po);
        std::vector<int> owner(64, -1);
        bool split = false;
        for (Id c = 0; c < m.cx.n_cells(); ++c)
        {
            if (owner[grp[c]] < 0) owner[grp[c]] = part[c];
            else if (owner[grp[c]] != part[c]) split = true;
        }
        const PartitionReport r = evaluate_partition(cg, part, np);
        std::printf("[grouped partition np=%d] 32 groups: split=%d imbalance %.3f cut faces %ld\n", np, int(split), r.imbalance, r.total_cut_faces);
        CHECK(!split);
        CHECK(r.imbalance < 1.0 + 1.0 * np / 32.0 + 0.05); // one group is 1/32 of the load
    }
    TEST_END();
}
