// Sub-complex extraction: closure, ownership, symmetric plans, and an in-process
// simulation of the edge reduce -> owner -> broadcast cycle used by the EMF assembly.
#include <algorithm>
#include <cstdio>
#include <set>

#include "TestUtil.h"
#include "uns/complex/Validate.h"
#include "uns/generate/HexMesh.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"

using namespace uns;

static int g_cells_of_edge(const CellComplex &cx, Id e)
{
    std::set<Id> s;
    for (auto k = cx.e2f.begin(e); k < cx.e2f.end(e); ++k)
        for (auto q = cx.f2c.begin(cx.e2f.idx[k]); q < cx.f2c.end(cx.e2f.idx[k]); ++q) s.insert(cx.f2c.idx[q]);
    return static_cast<int>(s.size());
}

static void run(const char *name, const Mesh &m, int np, PartMethod pm)
{
    Geometry geo = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, geo, CostModel());
    PartOptions po;
    po.method = pm;
    auto part = partition_cells(cg, np, po);
    Distributor dist(m, part, np);
    std::vector<LocalMesh> lm;
    for (int r = 0; r < np; ++r) lm.push_back(dist.build(r));

    long owned[4] = {0, 0, 0, 0};
    for (int r = 0; r < np; ++r)
    {
        const CellComplex &cx = lm[r].mesh.cx;
        TopologyReport t = validate_topology(cx);
        CHECK(t.d1d0_violations == 0 && t.d2d1_violations == 0 && t.bad_edges == 0 && t.bad_loops == 0);
        for (int d = 0; d < 4; ++d)
        {
            owned[d] += cx.ent[d].n_owned;
            for (Id i = 0; i < cx.ent[d].n_owned; ++i) CHECK(cx.ent[d].owner[i] == r);
            for (Id i = cx.ent[d].n_owned; i < cx.count(static_cast<Dim>(d)); ++i) CHECK(cx.ent[d].owner[i] != r);
        }
        // every edge of an owned cell sees all its incident cells locally
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
                for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke)
                {
                    const Id e = cx.f2e.idx[ke];
                    CHECK(g_cells_of_edge(cx, e) == g_cells_of_edge(m.cx, static_cast<Id>(cx.ent[1].gid[e])));
                }
    }
    for (int d = 0; d < 4; ++d) CHECK(owned[d] == m.cx.count(static_cast<Dim>(d)));

    // plans are symmetric: what r sends to q is what q expects from r
    for (int d = 0; d < 4; ++d)
        for (int r = 0; r < np; ++r)
            for (std::size_t i = 0; i < lm[r].plan.send[d].peers.size(); ++i)
            {
                const int q = lm[r].plan.send[d].peers[i];
                const auto &rq = lm[q].plan.recv[d];
                auto it = std::find(rq.peers.begin(), rq.peers.end(), r);
                CHECK(it != rq.peers.end());
                if (it == rq.peers.end()) continue;
                const std::size_t j = it - rq.peers.begin();
                CHECK(lm[r].plan.send[d].ptr[i + 1] - lm[r].plan.send[d].ptr[i] == rq.ptr[j + 1] - rq.ptr[j]);
                for (std::int64_t k = 0; k < rq.ptr[j + 1] - rq.ptr[j]; ++k)
                    CHECK(lm[r].mesh.cx.ent[d].gid[lm[r].plan.send[d].idx[lm[r].plan.send[d].ptr[i] + k]] ==
                          lm[q].mesh.cx.ent[d].gid[rq.idx[rq.ptr[j] + k]]);
            }

    // simulated edge assembly: each owned cell adds 1 to each of its edges
    std::vector<std::vector<double>> val(np);
    for (int r = 0; r < np; ++r)
    {
        const CellComplex &cx = lm[r].mesh.cx;
        val[r].assign(cx.n_edges(), 0.0);
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            std::set<Id> es;
            for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
                for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke) es.insert(cx.f2e.idx[ke]);
            for (Id e : es) val[r][e] += 1.0;
        }
    }
    for (int r = 0; r < np; ++r) // reduce alias -> owner
        for (std::size_t i = 0; i < lm[r].plan.reduce_send.peers.size(); ++i)
        {
            const int o = lm[r].plan.reduce_send.peers[i];
            const auto &rv = lm[o].plan.reduce_recv;
            const std::size_t j = std::find(rv.peers.begin(), rv.peers.end(), r) - rv.peers.begin();
            CHECK(j < rv.peers.size());
            for (std::int64_t k = lm[r].plan.reduce_send.ptr[i]; k < lm[r].plan.reduce_send.ptr[i + 1]; ++k)
                val[o][rv.idx[rv.ptr[j] + (k - lm[r].plan.reduce_send.ptr[i])]] += val[r][lm[r].plan.reduce_send.idx[k]];
        }
    for (int r = 0; r < np; ++r) // broadcast owner -> alias
        for (std::size_t i = 0; i < lm[r].plan.send[1].peers.size(); ++i)
        {
            const int q = lm[r].plan.send[1].peers[i];
            const auto &rq = lm[q].plan.recv[1];
            const std::size_t j = std::find(rq.peers.begin(), rq.peers.end(), r) - rq.peers.begin();
            for (std::int64_t k = lm[r].plan.send[1].ptr[i]; k < lm[r].plan.send[1].ptr[i + 1]; ++k)
                val[q][rq.idx[rq.ptr[j] + (k - lm[r].plan.send[1].ptr[i])]] = val[r][lm[r].plan.send[1].idx[k]];
        }
    long bad = 0;
    for (int r = 0; r < np; ++r)
    {
        const CellComplex &cx = lm[r].mesh.cx;
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
                for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke)
                {
                    const Id e = cx.f2e.idx[ke];
                    if (val[r][e] != g_cells_of_edge(m.cx, static_cast<Id>(cx.ent[1].gid[e]))) ++bad;
                }
    }
    CHECK(bad == 0);
    std::printf("[%s] np=%d %-7s ok (assembled-edge mismatches: %ld)\n", name, np, to_string(pm).c_str(), bad);
}

int main()
{
    BoxSpec b;
    b.nx = 7; b.ny = 6; b.nz = 5;
    Mesh box = build_mesh(make_box(b));
    CubedSphereSpec cs;
    cs.ncube = 5; cs.nr = 6; cs.n_solid = 2; cs.r_in = 0.8; cs.r_out = 6.0; cs.r_stretch = 1.2;
    Mesh shell = build_mesh(make_cubed_sphere_shell(cs));
    for (int np : {1, 2, 3, 5, 8})
        for (PartMethod pm : {PartMethod::Hilbert, PartMethod::Rcb})
        {
            run("box", box, np, pm);
            run("shell", shell, np, pm);
        }
    TEST_END();
}
