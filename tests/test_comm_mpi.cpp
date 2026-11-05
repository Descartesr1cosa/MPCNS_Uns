// MPI: edge assembly through the owner (reduce + broadcast) reproduces the global
// incident-cell count on every edge of an owned cell; the 2-cochain forward copy is exact.
#include <mpi.h>

#include <cstdio>
#include <set>

#include "uns/comm/Exchanger.h"
#include "uns/generate/HexMesh.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"

using namespace uns;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    CubedSphereSpec cs;
    cs.ncube = 5; cs.nr = 5; cs.n_solid = 1; cs.r_in = 0.8; cs.r_out = 6.0;
    Mesh m = build_mesh(make_cubed_sphere_shell(cs));
    Geometry geo = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, geo, CostModel());
    PartOptions po;
    po.method = PartMethod::Rcb;
    auto part = partition_cells(cg, np, po);
    Distributor dist(m, part, np);
    LocalMesh lm = dist.build(rank);
    const CellComplex &cx = lm.mesh.cx;

    // each owned cell adds its sign-free "1" to all its edges (distinct per cell)
    std::vector<double> E(cx.n_edges(), 0.0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        std::set<Id> es;
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
            for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke) es.insert(cx.f2e.idx[ke]);
        for (Id e : es) E[e] += 1.0;
    }
    comm::Exchanger ex(MPI_COMM_WORLD, lm);
    ex.assemble_edges(E.data());

    long bad = 0;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
            for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke)
            {
                const Id e = cx.f2e.idx[ke];
                std::set<Id> s;
                const Id ge = static_cast<Id>(cx.ent[1].gid[e]);
                for (auto k = m.cx.e2f.begin(ge); k < m.cx.e2f.end(ge); ++k)
                    for (auto q = m.cx.f2c.begin(m.cx.e2f.idx[k]); q < m.cx.f2c.end(m.cx.e2f.idx[k]); ++q) s.insert(m.cx.f2c.idx[q]);
                if (E[e] != static_cast<double>(s.size())) ++bad;
            }
    // forward copy of a face field: value = global id, ghosts must receive it
    std::vector<double> B(cx.n_faces());
    for (Id f = 0; f < cx.n_faces(); ++f) B[f] = f < cx.ent[2].n_owned ? static_cast<double>(cx.ent[2].gid[f]) : -1.0;
    ex.forward(Dim::Face, B.data());
    for (Id f = 0; f < cx.n_faces(); ++f) if (B[f] != static_cast<double>(cx.ent[2].gid[f])) ++bad;

    long total = 0;
    MPI_Allreduce(&bad, &total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::printf("np=%d mismatches=%ld  %s\n", np, total, total ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return total ? 1 : 0;
}
