// Mars branch verification (MPI, any number of ranks):
//  A  free stream far from the planet (Hall + ambipolar on, all faces inflow, hanging-node mesh): H+ and the induced flux stay at the level of
//     the physical sources (hydrogen corona collisions / photo-ionisation); div B stays at round-off
//  B  crustal dipole background: B_add = d1(A.dl) is divergence free and matches the analytic dipole field
//  C  point-implicit sources: stiff ion-neutral drag / thermal relaxation / gyration stay bounded with dt far above the explicit limits
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <fstream>

#include "TestUtil.h"
#include "diagnostics/DivB.h"
#include "driver/MarsSolver.h"
#include "uns/generate/HexMesh.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"
#include "uns/refine/Refine.h"

using namespace uns;
using namespace mhd;
using namespace mars;

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

Params base_params()
{
    Params p;
    for (const char *a : {"bc_default=inflow", "write_vtk=0", "write_fields=0", "log_every=0", "is_Hall_Efield=1", "is_Ambipolar_Efield=1",
                          "mars_lower_boundary_km=100"})
        p.set_override(a);
    return p;
}

Mesh box_mesh(double x0, double h, int n)
{
    BoxSpec bs;
    bs.nx = bs.ny = bs.nz = n;
    bs.lo = {x0, x0 * 0.3, -0.5 * n * h};
    bs.hi = {x0 + n * h, x0 * 0.3 + n * h, 0.5 * n * h};
    return build_mesh(make_box(bs));
}

// ------------------------------------------------------------------------------------------------ A
void test_freestream(MPI_Comm comm, int rank, int np)
{
    Mesh box = box_mesh(30.0, 0.5, 4);
    RefineSpec rs;
    rs.max_level = 1;
    Region reg;
    reg.kind = Region::Kind::Sphere;
    reg.a = {31.0, 10.0, 0.0, 0.9};
    reg.level = 1;
    rs.regions.push_back(reg);
    Mesh hang = build_refined_mesh(box, rs, std::vector<int>(box.cx.n_cells(), 0)).fine;
    LocalMesh lm = local_of(hang, np, rank);
    MarsSolver s(comm, lm, MarsParams::from_params(base_params()));
    s.initialize();
    std::vector<double> U0[NS];
    for (int k = 0; k < NS; ++k) U0[k] = s.U(k);
    const std::vector<double> Phi0 = s.phi();
    for (int i = 0; i < 6; ++i) s.step(1e9);
    const CellComplex &cx = lm.mesh.cx;
    double errH = 0, errT = 0, errB = 0;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        for (int k = 0; k < NS; ++k)
            for (int m = 0; m < 5; ++m)
            {
                const double d = std::abs(s.U(k)[5 * c + m] - U0[k][5 * c + m]);
                (k == MARS::H ? errH : errT) = std::max(k == MARS::H ? errH : errT, d / (k == MARS::H ? 1.0 : 1e-10));
            }
    for (Id f = 0; f < cx.ent[2].n_owned; ++f) errB = std::max(errB, std::abs(s.phi()[f] - Phi0[f]));
    double loc[3] = {errH, errT, errB}, g[3];
    MPI_Allreduce(loc, g, 3, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0) std::printf("[A freestream, np=%d] max change after 6 steps: H+ %.2e (abs), trace ions %.2e (rel to floor), induced flux %.2e\n", np, g[0], g[1], g[2]);
    CHECK(g[0] < 1e-4);   // observed ~4e-6: ion-neutral drag on the H corona (1.3e3 cm^-3) and photo-ionisation
    CHECK(std::isfinite(g[1]));
    CHECK(g[2] < 1e-7);
}

// ------------------------------------------------------------------------------------------------ B
void test_dipole(MPI_Comm comm, int rank, int np)
{
    const double L = 3.396e6, B0 = 3e-9, mu0 = 1.256637061e-6, m = 1.0e20;
    const char *file = "test_mars_dipole.txt";
    if (rank == 0)
    {
        std::ofstream f(file);
        f << "MARS_DIPOLES_V1 1\n0 0 0 0 0 " << m << " 0\n";
    }
    MPI_Barrier(comm);
    Params p = base_params();
    p.set_override("is_Crustal_Bfield=1");
    p.set_override(std::string("Crustal_Dipoles_File=") + file);
    p.set_override("Crustal_Dipoles_Cutoff_km=0");
    p.set_override("sw_Bx=0");
    p.set_override("sw_By=0");
    Mesh box = box_mesh(1.5, 0.1, 10);
    LocalMesh lm = local_of(box, np, rank);
    MarsSolver s(comm, lm, MarsParams::from_params(p));
    s.initialize();
    const CellComplex &cx = lm.mesh.cx;
    const Geometry geo = compute_geometry(lm.mesh);
    const double dv = div_b(cx, s.phi_add().data()).max_rel;
    double num = 0, den = 0;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        const Vec3 r = geo.cell_center[c] * L;
        const double rn = norm(r);
        const Vec3 rh = (1.0 / rn) * r;
        const Vec3 mv(0, 0, m);
        const Vec3 b = (mu0 / (4 * 3.14159265358979) / (rn * rn * rn) / B0) * (3.0 * dot(mv, rh) * rh - mv);
        const Vec3 bn(s.B_total_cell()[3 * c], s.B_total_cell()[3 * c + 1], s.B_total_cell()[3 * c + 2]);
        num += dot(bn - b, bn - b);
        den += dot(b, b);
    }
    double loc[3] = {dv, num, den}, g[3];
    MPI_Allreduce(loc, g, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(loc + 1, g + 1, 2, MPI_DOUBLE, MPI_SUM, comm);
    const double rel = std::sqrt(g[1] / g[2]);
    if (rank == 0) std::printf("[B dipole background, np=%d] max rel div B_add = %.2e, L2 error of the cell field vs the analytic dipole = %.2e\n", np, g[0], rel);
    CHECK(g[0] < 1e-12);
    CHECK(rel < 0.03);
    MPI_Barrier(comm);
    if (rank == 0) std::remove(file);
}

// ------------------------------------------------------------------------------------------------ C
// Point-implicit sources with huge steps (explicit limits are ~1e-6 .. 1e-3):
//  C1 dense low-altitude gas: nu dt >> 1, drag must damp O2+ to rest, pressure stays positive and finite
//  C2 far from the planet (no gas): gyration with Omega dt >> 1 conserves the mass momentum  sum rho_s u_s  and stays bounded
void test_stiff(MPI_Comm comm, int rank, int np)
{
    const CellComplex *cxp = nullptr;
    {
        Params p = base_params();
        p.set_override("is_Hall_Efield=0");
        p.set_override("bc_default=outflow");
        BoxSpec bs;
        bs.nx = bs.ny = bs.nz = 4;
        bs.lo = {1.03, -0.1, -0.1};
        bs.hi = {1.13, 0.1, 0.1};
        Mesh box = build_mesh(make_box(bs));
        LocalMesh lm = local_of(box, np, rank);
        cxp = &lm.mesh.cx;
        MarsSolver s(comm, lm, MarsParams::from_params(p));
        s.initialize();
        const int O2 = MARS::O2;
        for (Id c = 0; c < cxp->n_cells(); ++c)
        {
            double *u = &s.U_mutable(O2)[5 * c];
            u[0] = 1e3; u[1] = 1e3 * 0.5; u[2] = 0; u[3] = 0; u[4] = 1e-3 / 0.4 + 0.5 * 1e3 * 0.25;
        }
        s.apply_implicit_sources(1.0);
        const Geometry geo1 = compute_geometry(lm.mesh);
        double umax = 0, pmin = 1e300;
        bool finite = true;
        for (Id c = 0; c < cxp->ent[3].n_owned; ++c)
            for (int k = 0; k < NS; ++k)
            {
                const double *u = &s.U_mutable(k)[5 * c];
                finite = finite && std::isfinite(u[0]) && std::isfinite(u[1]) && std::isfinite(u[4]);
                if (k == O2 && norm(geo1.cell_center[c]) < 1.045)
                {
                    umax = std::max(umax, std::abs(u[1]) / u[0]);
                    pmin = std::min(pmin, 0.4 * (u[4] - 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0]));
                }
            }
        double loc[2] = {umax, -pmin}, g[2];
        MPI_Allreduce(loc, g, 2, MPI_DOUBLE, MPI_MAX, comm);
        if (rank == 0) std::printf("[C1 drag, np=%d] one implicit step dt = 1: max |u_x| of O2+ below 150 km %.3e (from 0.5), min p %.3e\n", np, g[0], -g[1]);
        CHECK(finite);
        CHECK(g[0] < 1e-2);
        CHECK(-g[1] > 0);
    }
    {
        Params p = base_params();
        p.set_override("is_Hall_Efield=0");
        p.set_override("bc_default=outflow");
        Mesh box = box_mesh(30.0, 0.5, 3);
        LocalMesh lm = local_of(box, np, rank);
        MarsSolver s(comm, lm, MarsParams::from_params(p));
        s.initialize();
        const CellComplex &cx = lm.mesh.cx;
        const double mass[3] = {MARS::ParticleMass[MARS::H], MARS::ParticleMass[MARS::O], MARS::ParticleMass[MARS::O2]};
        const int sp[3] = {MARS::H, MARS::O, MARS::O2};
        const double rho_ref = 1.0 / mass[0]; // number densities n_H : n_O : n_O2 = 1 : 0.5 : 0.25 (nondimensional)
        (void)rho_ref;
        for (Id c = 0; c < cx.n_cells(); ++c)
            for (int q = 0; q < 3; ++q)
            {
                const double n = 1.0 / double(1 << q);
                const double rho = n * mass[q] / mass[0];
                const double vel[3][3] = {{-15, 0, 0}, {2, 3, 1}, {-4, 1, 5}};
                double *u = &s.U_mutable(sp[q])[5 * c];
                u[0] = rho;
                for (int a = 0; a < 3; ++a) u[1 + a] = rho * vel[q][a];
                u[4] = 1.0 + 0.5 * rho * (vel[q][0] * vel[q][0] + vel[q][1] * vel[q][1] + vel[q][2] * vel[q][2]);
            }
        auto momentum = [&](double out[3]) {
            s.prepare();
            const CellComplex &cxx = lm.mesh.cx;
            out[0] = out[1] = out[2] = 0;
            for (Id c = 0; c < cxx.ent[3].n_owned; ++c)
                for (int k = 0; k < NS; ++k)
                {
                    for (int a = 0; a < 3; ++a) out[a] += s.U(k)[5 * c + 1 + a]; // mass momentum  sum_s rho_s u_s
                }
        };
        double m0[3], m1[3], g0[3], g1[3];
        momentum(m0);
        s.apply_implicit_sources(50.0); // Omega_H dt ~ 1e3
        momentum(m1);
        MPI_Allreduce(m0, g0, 3, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(m1, g1, 3, MPI_DOUBLE, MPI_SUM, comm);
        double vmax = 0;
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
            for (int k = 0; k < NS; ++k)
                for (int a = 0; a < 3; ++a) vmax = std::max(vmax, std::abs(s.U(k)[5 * c + 1 + a] / s.U(k)[5 * c]));
        double gv;
        MPI_Allreduce(&vmax, &gv, 1, MPI_DOUBLE, MPI_MAX, comm);
        double dm = 0, sc = 0;
        for (int a = 0; a < 3; ++a) { dm = std::max(dm, std::abs(g1[a] - g0[a])); sc = std::max(sc, std::abs(g0[a])); }
        if (rank == 0) std::printf("   sum rho u before (%.4e %.4e %.4e) after (%.4e %.4e %.4e)\n", g0[0], g0[1], g0[2], g1[0], g1[1], g1[2]);
        if (rank == 0) std::printf("[C2 gyration, np=%d] Omega dt >> 1: change of sum rho_s u_s = %.2e (relative), max |u| %.3f (initial max 15)\n", np, dm / sc, gv);
        CHECK(dm / sc < 3e-3);  // the only other effect is drag on the H corona (1.3e3 cm^-3): nu dt = 2.7e-5 * 50 = 1.3e-3
        CHECK(gv < 20.0);
    }
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    test_freestream(MPI_COMM_WORLD, rank, np);
    test_dipole(MPI_COMM_WORLD, rank, np);
    test_stiff(MPI_COMM_WORLD, rank, np);
    int fails = g_failures, all = 0;
    MPI_Allreduce(&fails, &all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", all ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return all ? 1 : 0;
}
