// Pseudo-time (steady) mode verification (MPI, any number of ranks):
//  A  faraday_update_edge_dt keeps d2 B = 0 for arbitrary positive per-edge steps (hanging-node mesh, arbitrary edge valence)
//  B  steady mode with all steps capped at the global minimum (lts_max_ratio = 1) reproduces the global-step first-order scheme
//  C  steady mode with a strongly varying step field: d2 B stays at round-off, the state stays finite and positive
#include <mpi.h>

#include <cmath>
#include <cstdio>

#include "TestUtil.h"
#include "diagnostics/DivB.h"
#include "driver/Solver.h"
#include "induction/Faraday.h"
#include "uns/generate/HexMesh.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"
#include "uns/refine/Refine.h"

using namespace uns;
using namespace mhd;
using namespace mercury;

namespace
{
LocalMesh local_of(const Mesh &m, int nranks, int rank)
{
    Geometry g = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, g, CostModel());
    PartOptions po;
    po.method = PartMethod::Rcb;
    Distributor d(m, partition_cells(cg, nranks, po), nranks);
    return d.build(rank);
}

Mesh hanging_box()
{
    BoxSpec bs;
    bs.nx = 6; bs.ny = 5; bs.nz = 5; bs.hi = {6, 5, 5};
    Mesh box = build_mesh(make_box(bs));
    RefineSpec rs; rs.max_level = 2;
    Region reg; reg.kind = Region::Kind::Sphere; reg.a = {2.2, 1.8, 2.1, 1.3}; reg.level = 2;
    rs.regions.push_back(reg);
    return build_refined_mesh(box, rs, std::vector<int>(box.cx.n_cells(), 0)).fine;
}

Config pseudo_config(bool steady, double ratio)
{
    Config c;
    c.sc = Scales::simple(5.0 / 3.0, 1.0);
    c.lorentz_force = true;
    c.time_order = 1;
    c.second_order = false;
    c.bc_default = "outflow";
    c.fsH = {1.0, 1.0}; c.fsNa = {0.05, 0.05};
    c.u_inf = {0.5, 0.1, -0.2};
    c.B_imf = {0.4, 0.3, 0.2};
    c.steady = steady;
    c.lts_max_ratio = ratio;
    c.write_vtk = false; c.write_fields = false; c.write_profile = false; c.log_every = 0;
    return c;
}

void fill(Solver &s)
{
    const CellComplex &cx = s.mesh().mesh.cx;
    const double g = 5.0 / 3.0;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        const Vec3 x = s.geometry().cell_center[c];
        const Prim h{1.0 + 0.3 * std::sin(0.9 * x.x) * std::cos(0.7 * x.y), 0.5 + 0.2 * std::sin(x.z), 0.1 * std::cos(x.x), -0.2 + 0.1 * std::sin(x.y), 1.0 + 0.2 * std::cos(x.x + x.z)};
        const Prim n{0.05 * (1 + 0.5 * std::sin(x.x)), 0.4, 0.1, -0.2, 0.05};
        to_cons(h, g, &s.state().U[0][5 * static_cast<std::size_t>(c)]);
        to_cons(n, g, &s.state().U[1][5 * static_cast<std::size_t>(c)]);
    }
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    const Mesh m = hanging_box();
    LocalMesh lm = local_of(m, np, rank);
    const CellComplex &cx = lm.mesh.cx;

    // ------------------------------------------------------------------ A
    {
        std::vector<double> Phi(cx.n_faces(), 0.0), E(cx.n_edges()), dt(cx.n_edges());
        for (Id e = 0; e < cx.n_edges(); ++e) E[e] = std::sin(0.37 * static_cast<double>(cx.ent[1].gid[e]));
        // alias edges carry the owner's sign: the line integral of a cochain flips with the orientation, so use the sign-corrected value
        for (Id e = 0; e < cx.n_edges(); ++e) E[e] *= cx.ent[1].alias_sign[e];
        for (Id e = 0; e < cx.n_edges(); ++e) dt[e] = 0.2 + 3.0 * std::abs(std::cos(1.3 * static_cast<double>(cx.ent[1].gid[e])));
        const auto active = active_faces(cx);
        // start from a divergence-free field d1 A
        std::vector<double> A(cx.n_edges());
        for (Id e = 0; e < cx.n_edges(); ++e) A[e] = cx.ent[1].alias_sign[e] * std::cos(0.11 * static_cast<double>(cx.ent[1].gid[e]));
        for (Id f = 0; f < cx.n_faces(); ++f)
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) Phi[f] += cx.f2e.sgn[k] * A[cx.f2e.idx[k]];
        double worst = 0;
        for (int it = 0; it < 20; ++it)
        {
            faraday_update_edge_dt(cx, active, E.data(), dt.data(), Phi.data());
            worst = std::max(worst, mhd::div_b(cx, Phi.data()).max_rel);
        }
        double g = 0;
        MPI_Allreduce(&worst, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[A edge-step CT, np=%d] max relative d2 B over 20 updates = %.2e\n", np, g);
        CHECK(g < 1e-12);
    }

    // ------------------------------------------------------------------ B
    {
        Solver a(MPI_COMM_WORLD, lm, pseudo_config(false, 0.0)), b(MPI_COMM_WORLD, lm, pseudo_config(true, 1.0));
        a.initialize(); b.initialize();
        fill(a); fill(b);
        for (int i = 0; i < 8; ++i) { a.step(1e9); b.step(1e9); }
        double err = 0, ref = 0;
        for (int k = 0; k < 2; ++k)
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
                for (int q = 0; q < 5; ++q)
                {
                    const std::size_t i = 5 * static_cast<std::size_t>(c) + q;
                    err = std::max(err, std::abs(a.state().U[k][i] - b.state().U[k][i]));
                    ref = std::max(ref, std::abs(a.state().U[k][i]));
                }
        for (Id f = 0; f < cx.ent[2].n_owned; ++f) err = std::max(err, std::abs(a.state().Phi[f] - b.state().Phi[f]));
        double w[2] = {err, ref}, g[2];
        MPI_Allreduce(w, g, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[B ratio-1 pseudo time vs global step, np=%d] max difference after 8 steps = %.2e (|U| up to %.2f)\n", np, g[0], g[1]);
        CHECK(g[0] < 1e-12);
    }

    // ------------------------------------------------------------------ C
    {
        Solver s(MPI_COMM_WORLD, lm, pseudo_config(true, 0.0));
        s.initialize();
        fill(s);
        double mn = 1e300;
        for (int i = 0; i < 40; ++i) s.step(1e9);
        const Diagnostics d = s.diagnostics();
        mn = d.min_rho;
        if (rank == 0) std::printf("[C local steps, np=%d] 40 iterations: div B = %.2e, min rho = %.3f, step range [%.2e, %.2e]\n", np, d.div_b, mn, d.dt_min, d.dt_max);
        CHECK(d.div_b < 1e-12);
        CHECK(std::isfinite(mn) && mn > 0.05);
        CHECK(d.dt_max > 1.2 * d.dt_min); // the step field really varies
    }
    int fail = g_failures, gfail = 0;
    MPI_Allreduce(&fail, &gfail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", gfail ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return gfail ? 1 : 0;
}
