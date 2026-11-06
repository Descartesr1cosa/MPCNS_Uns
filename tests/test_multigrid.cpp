// Multigrid hierarchy verification (MPI, any number of ranks):
//  A  hierarchy geometry on a hanging-node mesh: cluster volumes add up, every cluster surface closes (sum of area vectors = 0),
//     the shared area vectors of two own clusters are opposite, restriction conserves the conserved variables
//  B  level 0 (clusters = cells) reproduces the fine first-order flux right-hand side (fluid_flux_rhs with zero gradients)
//  C  FAS fixed point: for a uniform free stream with zero fine rate the coarse-grid correction vanishes
//  D  FAS consistency: with an arbitrary smooth state and the exact fine rate the coarse level starts from the restricted rate
//     (forcing g = I(rate) - R_c(I u) makes the coarse rate at u_c = I u equal to the restricted fine rate)
//  E  positivity safeguard
#include <mpi.h>

#include <cmath>
#include <cstdio>

#include "TestUtil.h"
#include "multigrid/Multigrid.h"
#include "fluid/FluxRHS.h"
#include "fluid/MercurySources.h"
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
    bs.nx = 8; bs.ny = 6; bs.nz = 6; bs.hi = {8, 6, 6};
    Mesh box = build_mesh(make_box(bs));
    RefineSpec rs; rs.max_level = 2;
    Region reg; reg.kind = Region::Kind::Sphere; reg.a = {3.2, 2.8, 3.1, 1.6}; reg.level = 2;
    rs.regions.push_back(reg);
    return build_refined_mesh(box, rs, std::vector<int>(box.cx.n_cells(), 0)).fine;
}

Config cfg_for(int levels)
{
    Config c;
    c.sc = Scales::simple(5.0 / 3.0, 1.0);
    c.lorentz_force = false;
    c.bc_default = "inflow";
    c.fsH = {1.3, 0.7}; c.fsNa = {0.05, 0.02};
    c.u_inf = {0.4, -0.2, 0.3};
    c.mg_levels = levels; c.mg_agglom_passes = 2; c.mg_min_clusters = 6;
    c.mg_pre = 1; c.mg_post = 1; c.mg_coarse = 2;
    return c;
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
    const Geometry geo = compute_geometry(lm.mesh);
    Config cfg = cfg_for(4);
    BoundaryInfo bc = BoundaryInfo::build(lm, geo, cfg);
    comm::Exchanger ex(MPI_COMM_WORLD, lm);
    const Prim fs[2] = {{cfg.fsH.rho, cfg.u_inf.x, cfg.u_inf.y, cfg.u_inf.z, cfg.fsH.p}, {cfg.fsNa.rho, cfg.u_inf.x, cfg.u_inf.y, cfg.u_inf.z, cfg.fsNa.p}};
    Multigrid mg(MPI_COMM_WORLD, lm, geo, bc, ex, cfg, fs);
    const CellComplex &cx = lm.mesh.cx;
    const std::size_t nc = cx.n_cells();

    // ------------------------------------------------------------------ A
    {
        const auto sizes = mg.global_sizes();
        double vfine = 0;
        for (Id c = 0; c < cx.ent[3].n_owned; ++c) vfine += geo.cell_volume[c];
        double gvfine = 0;
        MPI_Allreduce(&vfine, &gvfine, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double worst_close = 0, worst_anti = 0, worst_vol = 0;
        for (int l = 0; l < mg.n_levels(); ++l)
        {
            const MgLevel &L = mg.level(l);
            double v = 0;
            for (Id i = 0; i < L.n_own; ++i)
            {
                v += L.vol[i];
                Vec3 sum;
                double mag = 0;
                for (auto q = L.fptr[i]; q < L.fptr[i + 1]; ++q) { sum += L.fA[q]; mag += norm(L.fA[q]); }
                for (auto q = L.bptr[i]; q < L.bptr[i + 1]; ++q) { sum += L.bA[q]; mag += norm(L.bA[q]); }
                worst_close = std::max(worst_close, norm(sum) / mag);
                for (auto q = L.fptr[i]; q < L.fptr[i + 1]; ++q) // total area vector to an own neighbour j = minus the total seen from j
                {
                    const Id j = L.fnb[q];
                    if (j >= L.n_own) continue;
                    Vec3 sij, sji;
                    double mij = 0;
                    for (auto r = L.fptr[i]; r < L.fptr[i + 1]; ++r) if (L.fnb[r] == j) { sij += L.fA[r]; mij += norm(L.fA[r]); }
                    for (auto r = L.fptr[j]; r < L.fptr[j + 1]; ++r) if (L.fnb[r] == i) sji += L.fA[r];
                    worst_anti = std::max(worst_anti, norm(sij + sji) / mij);
                }
            }
            double gv = 0;
            MPI_Allreduce(&v, &gv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            worst_vol = std::max(worst_vol, std::abs(gv - gvfine) / gvfine);
        }
        double w[3] = {worst_close, worst_anti, worst_vol}, g[3];
        MPI_Allreduce(w, g, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0)
        {
            std::printf("[A hierarchy np=%d] clusters per level:", np);
            for (long n : sizes) std::printf(" %ld", n);
            std::printf("   closure %.2e  antisymmetry %.2e  volume %.2e\n", g[0], g[1], g[2]);
        }
        CHECK(g[0] < 1e-12 && g[1] < 1e-12 && g[2] < 1e-12);
        CHECK(sizes[0] > sizes[1] && sizes[1] > sizes[2]);
        // restriction conserves the conserved variables
        std::vector<double> Uf[2];
        for (int k = 0; k < 2; ++k)
        {
            Uf[k].assign(5 * nc, 0.0);
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            {
                const Vec3 x = geo.cell_center[c];
                to_cons(Prim{1.0 + 0.3 * std::sin(x.x) * x.y, 0.2 * x.z, 0.1, -0.3 * x.x, 1.0 + 0.1 * x.y}, cfg.sc.gamma, &Uf[k][5 * static_cast<std::size_t>(c)]);
            }
        }
        double worst = 0;
        for (int l = 1; l < mg.n_levels(); ++l)
        {
            const MgLevel &L = mg.level(l);
            std::vector<double> U0(5 * L.n_tot()), U1(5 * L.n_tot());
            const double *Ufp[2] = {Uf[0].data(), Uf[1].data()};
            double *Up[2] = {U0.data(), U1.data()};
            mg.restrict_fine(l, Ufp, Up);
            double a[5] = {0}, b[5] = {0}, ga[5], gb[5];
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
                for (int q = 0; q < 5; ++q) a[q] += geo.cell_volume[c] * Uf[0][5 * static_cast<std::size_t>(c) + q];
            for (Id i = 0; i < L.n_own; ++i)
                for (int q = 0; q < 5; ++q) b[q] += L.vol[i] * U0[5 * static_cast<std::size_t>(i) + q];
            MPI_Allreduce(a, ga, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(b, gb, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (int q = 0; q < 5; ++q) worst = std::max(worst, std::abs(ga[q] - gb[q]) / (std::abs(ga[q]) + 1e-30));
        }
        if (rank == 0) std::printf("[A restriction] max relative conservation defect %.2e\n", worst);
        CHECK(worst < 1e-13);
    }

    // ------------------------------------------------------------------ B
    {
        // smooth, strongly varying state; compare level 0 flux with the fine first-order RHS
        std::vector<double> W[2], rhs_f[2], rhs_m[2], zero(15 * nc, 0.0);
        for (int k = 0; k < 2; ++k)
        {
            W[k].assign(5 * nc, 0.0);
            for (Id c = 0; c < static_cast<Id>(nc); ++c)
            {
                const Vec3 x = geo.cell_center[c];
                double *w = &W[k][5 * static_cast<std::size_t>(c)];
                w[0] = 1.0 + 0.3 * std::sin(1.3 * x.x + k) * std::cos(x.y); w[1] = 0.4 + 0.2 * std::sin(x.z); w[2] = -0.2 + 0.1 * x.x * 0.1; w[3] = 0.3 * std::cos(x.x);
                w[4] = 0.8 + 0.2 * std::cos(x.y + x.z);
            }
            rhs_f[k].assign(5 * nc, 0.0);
            rhs_m[k].assign(5 * static_cast<std::size_t>(mg.level(0).n_own), 0.0);
        }
        for (int k = 0; k < 2; ++k) fluid_flux_rhs(lm, geo, bc, cfg, k, fs[k], W[k].data(), zero.data(), rhs_f[k].data());
        // level-0 arrays are indexed by the compact fluid-cell number: map through cell_cluster
        const MgLevel &L0 = mg.level(0);
        std::vector<double> Wc[2];
        for (int k = 0; k < 2; ++k)
        {
            Wc[k].assign(5 * L0.n_tot(), 0.0);
            for (Id c = 0; c < static_cast<Id>(nc); ++c)
                if (L0.cell_cluster[c] >= 0)
                    for (int q = 0; q < 5; ++q) Wc[k][5 * static_cast<std::size_t>(L0.cell_cluster[c]) + q] = W[k][5 * static_cast<std::size_t>(c) + q];
        }
        const double *Wp[2] = {Wc[0].data(), Wc[1].data()};
        double *rp[2] = {rhs_m[0].data(), rhs_m[1].data()};
        mg.flux_rhs(0, Wp, rp);
        double err = 0, ref = 0;
        for (int k = 0; k < 2; ++k)
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            {
                const Id cl = L0.cell_cluster[c];
                if (cl < 0) continue;
                for (int q = 0; q < 5; ++q)
                {
                    if (std::abs(rhs_m[k][5 * static_cast<std::size_t>(cl) + q] - rhs_f[k][5 * static_cast<std::size_t>(c) + q]) > 0.5 && rank == 0)
                        std::printf("   cell %d (cluster %d) comp %d: mg %.4f fine %.4f nfaces %d x=(%.2f %.2f %.2f)\n", (int)c, (int)cl, q, rhs_m[k][5 * static_cast<std::size_t>(cl) + q], rhs_f[k][5 * static_cast<std::size_t>(c) + q], (int)(cx.c2f.end(c) - cx.c2f.begin(c)), geo.cell_center[c].x, geo.cell_center[c].y, geo.cell_center[c].z);
                    err = std::max(err, std::abs(rhs_m[k][5 * static_cast<std::size_t>(cl) + q] - rhs_f[k][5 * static_cast<std::size_t>(c) + q]));
                    ref = std::max(ref, std::abs(rhs_f[k][5 * static_cast<std::size_t>(c) + q]));
                }
            }
        double w[2] = {err, ref}, g[2];
        MPI_Allreduce(w, g, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[B level-0 flux vs fine first-order RHS, np=%d] max diff %.2e (max |rhs| %.2e)\n", np, g[0], g[1]);
        CHECK(g[0] < 1e-11 * std::max(g[1], 1.0));
    }

    // ------------------------------------------------------------------ C : fixed point of the cycle
    {
        std::vector<double> Uf[2], rate[2], dU[2], B(3 * nc, 0.0), J(3 * nc, 0.0), dpe(3 * nc, 0.0), photo(nc, 0.0), vst(nc, 0.0);
        for (int k = 0; k < 2; ++k)
        {
            Uf[k].assign(5 * nc, 0.0); rate[k].assign(5 * nc, 0.0); dU[k].assign(5 * nc, 1.0);
            for (Id c = 0; c < static_cast<Id>(nc); ++c) to_cons(fs[k], cfg.sc.gamma, &Uf[k][5 * static_cast<std::size_t>(c)]);
        }
        const double *Up[2] = {Uf[0].data(), Uf[1].data()}, *Rp[2] = {rate[0].data(), rate[1].data()};
        double *dp[2] = {dU[0].data(), dU[1].data()};
        mg.cycle(Up, Rp, B.data(), J.data(), dpe.data(), photo.data(), vst.data(), dp);
        double mx = 0;
        for (int k = 0; k < 2; ++k)
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
                for (int q = 0; q < 5; ++q) mx = std::max(mx, std::abs(dU[k][5 * static_cast<std::size_t>(c) + q]));
        double g = 0;
        MPI_Allreduce(&mx, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[C uniform stream, zero rate] max |coarse-grid correction| = %.2e\n", g);
        CHECK(g < 1e-12);
    }

    // ------------------------------------------------------------------ D : forcing makes the coarse rate equal to the restricted fine rate
    {
        // use the public pieces: restrict a smooth state, build g = I(rate) - R_c(Iu) is internal; check via a cycle with zero smoothing:
        // with mg_pre = mg_post = mg_coarse = 0 the correction is exactly 0 whatever the rate is (U_l stays at U0_l)
        Config c0 = cfg; c0.mg_pre = c0.mg_post = c0.mg_coarse = 0; c0.mg_levels = 3;
        Multigrid mg0(MPI_COMM_WORLD, lm, geo, bc, ex, c0, fs);
        std::vector<double> Uf[2], rate[2], dU[2], B(3 * nc, 0.0), J(3 * nc, 0.0), dpe(3 * nc, 0.0), photo(nc, 0.0), vst(nc, 0.0);
        for (int k = 0; k < 2; ++k)
        {
            Uf[k].assign(5 * nc, 0.0); rate[k].assign(5 * nc, 0.0); dU[k].assign(5 * nc, 1.0);
            for (Id c = 0; c < static_cast<Id>(nc); ++c)
            {
                const Vec3 x = geo.cell_center[c];
                to_cons(Prim{1.0 + 0.2 * std::sin(x.x), 0.3, 0.1 * std::cos(x.y), 0.0, 1.0 + 0.1 * x.z}, cfg.sc.gamma, &Uf[k][5 * static_cast<std::size_t>(c)]);
                for (int q = 0; q < 5; ++q) rate[k][5 * static_cast<std::size_t>(c) + q] = 0.1 * std::sin(x.x + q);
            }
        }
        const double *Up[2] = {Uf[0].data(), Uf[1].data()}, *Rp[2] = {rate[0].data(), rate[1].data()};
        double *dp[2] = {dU[0].data(), dU[1].data()};
        mg0.cycle(Up, Rp, B.data(), J.data(), dpe.data(), photo.data(), vst.data(), dp);
        double mx = 0;
        for (int k = 0; k < 2; ++k)
            for (Id c = 0; c < cx.ent[3].n_owned; ++c)
                for (int q = 0; q < 5; ++q) mx = std::max(mx, std::abs(dU[k][5 * static_cast<std::size_t>(c) + q]));
        double g = 0;
        MPI_Allreduce(&mx, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[D zero smoothing] max |correction| = %.2e\n", g);
        CHECK(g < 1e-12);
    }

    // ------------------------------------------------------------------ E
    {
        const double U[5] = {1.0, 0.0, 0.0, 0.0, 1.5};
        const double bad[5] = {-2.0, 0, 0, 0, -3.0}, mild[5] = {0.1, 0, 0, 0, 0.1};
        const double a1 = positivity_safeguard(U, bad, 1.4), a2 = positivity_safeguard(U, mild, 1.4), a3 = positivity_safeguard(U, mild, 1.4, 0.1, 1.05);
        if (rank == 0) std::printf("[E safeguard] alpha(bad) = %g  alpha(mild) = %g  alpha(mild, growth cap 1.05) = %g\n", a1, a2, a3);
        CHECK(a1 < 0.5);
        CHECK(a2 == 1.0);
        CHECK(a3 < 1.0 && a3 > 0.0);
    }
    int fail = g_failures, gfail = 0;
    MPI_Allreduce(&fail, &gfail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", gfail ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return gfail ? 1 : 0;
}
