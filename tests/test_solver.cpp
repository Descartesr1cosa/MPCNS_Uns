// Solver verification (MPI, any number of ranks):
//  A  free-stream preservation on a hanging-node polyhedral mesh (all faces inflow)       -> round-off
//  B  smooth density wave: observed order of the second-order scheme vs first order
//  C  Sod shock tube vs the exact Riemann solution
//  D  Alfven wave packet (Lorentz force + constrained transport): speed, amplitude, sign
//  E  Mercury-like shell (dipole, wall, far field, solid core): np-independence and div B
#include <mpi.h>

#include <cmath>
#include <cstdio>

#include "TestUtil.h"
#include "driver/Solver.h"
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

Config box_config(double gamma, double inver_MA2)
{
    Config c;
    c.sc = Scales::simple(gamma, inver_MA2);
    c.lorentz_force = inver_MA2 > 0;
    c.write_vtk = false;
    c.write_fields = false;
    c.log_every = 0;
    c.bc_default = "outflow";
    return c;
}

Mesh tube(int nx, double x0, double x1, int ncross)
{
    const double h = (x1 - x0) / nx;
    BoxSpec b;
    b.nx = nx; b.ny = ncross; b.nz = ncross;
    b.lo = {x0, 0, 0}; b.hi = {x1, ncross * h, ncross * h};
    return build_mesh(make_box(b));
}

// set a primitive state field from a function of the cell centre (all local cells, incl. ghosts)
template <class F> void set_state(Solver &s, int species, double gamma, F &&fn)
{
    const CellComplex &cx = s.mesh().mesh.cx;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        const Prim w = fn(s.geometry().cell_center[c]);
        to_cons(w, gamma, &s.state().U[species][5 * static_cast<std::size_t>(c)]);
    }
}

// ---------------------------------------------------------------------------------------------- A
void test_freestream(MPI_Comm comm, int rank, int np)
{
    BoxSpec bs;
    bs.nx = 5; bs.ny = 4; bs.nz = 4; bs.hi = {5, 4, 4};
    Mesh box = build_mesh(make_box(bs));
    RefineSpec rs; rs.max_level = 2;
    Region reg; reg.kind = Region::Kind::Sphere; reg.a = {2.2, 1.8, 2.1, 1.3}; reg.level = 2;
    rs.regions.push_back(reg);
    Mesh hang = build_refined_mesh(box, rs, std::vector<int>(box.cx.n_cells(), 0)).fine;
    Config c = box_config(5.0 / 3.0, 1.0);
    c.bc_default = "inflow";
    c.fsH = {1.3, 0.7}; c.fsNa = {0.01, 0.005};
    c.u_inf = {0.4, -0.2, 0.3};
    c.B_imf = {0.3, 0.5, -0.2};
    LocalMesh lm = local_of(hang, np, rank);
    Solver s(comm, lm, c);
    s.initialize();
    const State s0 = s.state();
    for (int i = 0; i < 6; ++i) s.step(1e9);
    const CellComplex &cx = lm.mesh.cx;
    double err = 0;
    for (Id cc = 0; cc < cx.ent[3].n_owned; ++cc)
        for (int k = 0; k < 2; ++k)
            for (int m = 0; m < 5; ++m) err = std::max(err, std::abs(s.state().U[k][5 * cc + m] - s0.U[k][5 * cc + m]));
    for (Id f = 0; f < cx.n_faces(); ++f)
        if (f < cx.ent[2].n_owned) err = std::max(err, std::abs(s.state().Phi[f] - s0.Phi[f]));
    double g = 0;
    MPI_Allreduce(&err, &g, 1, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0) std::printf("[A freestream, hanging mesh, np=%d] max change after 6 steps = %.2e\n", np, g);
    CHECK(g < 1e-12);
}

// ---------------------------------------------------------------------------------------------- B
double wave_error(MPI_Comm comm, int rank, int np, int nx, bool second_order)
{
    const double gamma = 1.4, sigma = 0.07, x0 = 0.35, amp = 0.2, uadv = 1.0, T = 0.25;
    Config c = box_config(gamma, 0.0);
    c.second_order = second_order;
    c.cfl = 0.6;
    c.t_end = T;
    c.bc_by_tag["xmin"] = "inflow";
    c.fsH = {1.0, 1.0}; c.fsNa = {1e-6, 1e-6}; c.u_inf = {uadv, 0, 0};
    Mesh m = tube(nx, 0.0, 1.0, 2);
    LocalMesh lm = local_of(m, np, rank);
    Solver s(comm, lm, c);
    s.initialize();
    auto rho0 = [&](double x) { return 1.0 + amp * std::exp(-(x - x0) * (x - x0) / (2 * sigma * sigma)); };
    set_state(s, 0, gamma, [&](const Vec3 &x) { return Prim{rho0(x.x), uadv, 0, 0, 1.0}; });
    while (s.time() < T * (1 - 1e-12)) s.step(T - s.time());
    s.diagnostics(); // refreshes primitives
    double e2 = 0, vol = 0;
    const CellComplex &cx = lm.mesh.cx;
    for (Id cc = 0; cc < cx.ent[3].n_owned; ++cc)
    {
        const double d = s.prim(0)[5 * cc] - rho0(s.geometry().cell_center[cc].x - uadv * T);
        e2 += d * d * s.geometry().cell_volume[cc];
        vol += s.geometry().cell_volume[cc];
    }
    double loc[2] = {e2, vol}, g[2];
    MPI_Allreduce(loc, g, 2, MPI_DOUBLE, MPI_SUM, comm);
    return std::sqrt(g[0] / g[1]);
}

void test_wave(MPI_Comm comm, int rank, int np)
{
    double e2[3], e1[3];
    const int ns[3] = {40, 80, 160};
    for (int i = 0; i < 3; ++i)
    {
        e2[i] = wave_error(comm, rank, np, ns[i], true);
        e1[i] = wave_error(comm, rank, np, ns[i], false);
    }
    if (rank == 0)
    {
        std::printf("[B density wave] L2 error  2nd order: %.3e %.3e %.3e (orders %.2f %.2f)   1st order: %.3e %.3e %.3e (orders %.2f %.2f)\n", e2[0],
                    e2[1], e2[2], std::log2(e2[0] / e2[1]), std::log2(e2[1] / e2[2]), e1[0], e1[1], e1[2], std::log2(e1[0] / e1[1]),
                    std::log2(e1[1] / e1[2]));
    }
    CHECK(std::log2(e2[1] / e2[2]) > 1.5);
    CHECK(std::log2(e1[1] / e1[2]) < 1.2);
    CHECK(e2[2] < 0.3 * e1[2]);
}

// ---------------------------------------------------------------------------------------------- C
double sod_exact_rho(double x, double t)
{
    const double g = 1.4, rL = 1.0, pL = 1.0, rR = 0.125, ps = 0.30313, us = 0.92745, rsL = 0.42632, rsR = 0.26557, S = 1.75216;
    const double cL = std::sqrt(g * pL / rL), x0 = 0.5, csL = std::sqrt(g * ps / rsL);
    const double xh = x0 - cL * t, xt = x0 + (us - csL) * t, xc = x0 + us * t, xs = x0 + S * t;
    if (x < xh) return rL;
    if (x < xt)
    {
        const double u = 2.0 / (g + 1) * (cL + (x - x0) / t), c = cL - 0.5 * (g - 1) * u;
        return rL * std::pow(c / cL, 2.0 / (g - 1));
    }
    if (x < xc) return rsL;
    if (x < xs) return rsR;
    return rR;
}

void test_sod(MPI_Comm comm, int rank, int np)
{
    const double gamma = 1.4, T = 0.2;
    Config c = box_config(gamma, 0.0);
    c.t_end = T;
    c.fsH = {1.0, 1.0}; c.fsNa = {1e-6, 1e-6}; c.u_inf = {0, 0, 0};
    Mesh m = tube(200, 0.0, 1.0, 2);
    LocalMesh lm = local_of(m, np, rank);
    Solver s(comm, lm, c);
    s.initialize();
    set_state(s, 0, gamma, [&](const Vec3 &x) { return x.x < 0.5 ? Prim{1.0, 0, 0, 0, 1.0} : Prim{0.125, 0, 0, 0, 0.1}; });
    while (s.time() < T * (1 - 1e-12)) s.step(T - s.time());
    s.diagnostics();
    double l1 = 0, vol = 0, mx = 0;
    const CellComplex &cx = lm.mesh.cx;
    for (Id cc = 0; cc < cx.ent[3].n_owned; ++cc)
    {
        const double d = std::abs(s.prim(0)[5 * cc] - sod_exact_rho(s.geometry().cell_center[cc].x, T));
        l1 += d * s.geometry().cell_volume[cc];
        vol += s.geometry().cell_volume[cc];
        mx = std::max(mx, s.prim(0)[5 * cc]);
    }
    double loc[2] = {l1, vol}, g[2], gmx;
    MPI_Allreduce(loc, g, 2, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&mx, &gmx, 1, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0) std::printf("[C Sod, 200 cells] L1 density error = %.3e, max rho = %.4f (exact 1)\n", g[0] / g[1], gmx);
    CHECK(g[0] / g[1] < 0.012);
    CHECK(gmx < 1.01);
}

// ---------------------------------------------------------------------------------------------- D
void test_alfven(MPI_Comm comm, int rank, int np)
{
    const double eps = 0.05, sigma = 0.4, T = 3.0;
    Config c = box_config(5.0 / 3.0, 1.0);
    c.t_end = T;
    c.fsH = {1.0, 0.1}; c.fsNa = {1e-6, 1e-7}; c.u_inf = {0, 0, 0};
    c.B_imf = {1.0, 0, 0};
    Mesh m = tube(240, -6.0, 6.0, 2);
    LocalMesh lm = local_of(m, np, rank);
    Solver s(comm, lm, c);
    s.initialize();
    set_state(s, 0, c.sc.gamma, [&](const Vec3 &x) { return Prim{1.0, 0, eps * std::exp(-x.x * x.x / (2 * sigma * sigma)), 0, 0.1}; });
    while (s.time() < T * (1 - 1e-12)) s.step(T - s.time());
    s.diagnostics();
    const CellComplex &cx = lm.mesh.cx;
    struct { double v; double x; } uy{-1e300, 0}, by{1e300, 0};
    for (Id cc = 0; cc < cx.ent[3].n_owned; ++cc)
    {
        const double x = s.geometry().cell_center[cc].x;
        if (x < 0.5) continue; // right-going packet only
        if (s.prim(0)[5 * cc + 2] > uy.v) uy = {s.prim(0)[5 * cc + 2], x};
        if (s.B_induced_cell()[3 * cc + 1] < by.v) by = {s.B_induced_cell()[3 * cc + 1], x};
    }
    // global max of u_y and min of B_y with their positions
    struct P { double v; int r; } in{uy.v, rank}, out;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    double xu = uy.x, xb = by.x, gxu, gxb;
    MPI_Bcast(&xu, 1, MPI_DOUBLE, out.r, comm);
    gxu = xu;
    P inb{-by.v, rank}, outb;
    MPI_Allreduce(&inb, &outb, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
    MPI_Bcast(&xb, 1, MPI_DOUBLE, outb.r, comm);
    gxb = xb;
    const double vA = 1.0;
    if (rank == 0)
        std::printf("[D Alfven packet, np=%d] right-going: u_y peak %.4f at x=%.3f (expect %.4f at %.3f);  B_y trough %.4f at x=%.3f (expect %.4f)\n", np,
                    out.v, gxu, eps / 2, vA * T, -outb.v, gxb, -eps / 2);
    CHECK(std::abs(out.v / (eps / 2) - 1.0) < 0.08);
    CHECK(std::abs(gxu - vA * T) < 0.12);
    CHECK(std::abs(-outb.v / (-eps / 2) - 1.0) < 0.08); // B_y = -u_y for the right-going wave with B0 along +x
    CHECK(std::abs(gxb - vA * T) < 0.12);
    const Diagnostics d = s.diagnostics();
    CHECK(d.div_b < 1e-13);
}

// ---------------------------------------------------------------------------------------------- E
void test_shell(MPI_Comm comm, int rank, int np)
{
    CubedSphereSpec cs;
    cs.ncube = 5; cs.nr = 7; cs.n_solid = 2; cs.r_in = 0.7; cs.r_out = 10; cs.r_stretch = 1.15;
    Mesh m = build_mesh(make_cubed_sphere_shell(cs));
    Params p;
    Config c = Config::from_params(p);
    c.sc = Scales::from_params(p);
    c.fsH = {1.0, 0.0089}; c.fsNa = {1e-6, 1e-8};
    c.u_inf = {-1.0, 0.12, 0};
    c.B_imf = {-15.2, 8.4, -8.5};
    c.dipole = true; c.dipole_m = {0, 0, -195}; c.dipole_c = {0, 0, 0.2};
    c.write_vtk = false; c.write_fields = false; c.log_every = 0;
    const int nsteps = 8;

    const LocalMesh lm_ref = local_of(m, 1, 0); // the Solver keeps a reference to its mesh
    Solver ref(MPI_COMM_SELF, lm_ref, c);
    ref.initialize();
    for (int i = 0; i < nsteps; ++i) ref.step(1e9);
    const Diagnostics dref = ref.diagnostics();

    LocalMesh lm = local_of(m, np, rank);
    Solver s(comm, lm, c);
    s.initialize();
    for (int i = 0; i < nsteps; ++i) s.step(1e9);
    const Diagnostics d = s.diagnostics();
    const CellComplex &cx = lm.mesh.cx;
    double err = 0, scale = 0;
    for (Id cc = 0; cc < cx.ent[3].n_owned; ++cc)
        for (int k = 0; k < 2; ++k)
            for (int q = 0; q < 5; ++q)
            {
                const double a = s.state().U[k][5 * cc + q], b = ref.state().U[k][5 * static_cast<std::size_t>(cx.ent[3].gid[cc]) + q];
                err = std::max(err, std::abs(a - b));
                scale = std::max(scale, std::abs(b));
            }
    double perr = 0;
    for (Id f = 0; f < cx.ent[2].n_owned; ++f) perr = std::max(perr, std::abs(s.state().Phi[f] - ref.state().Phi[cx.ent[2].gid[f]]));
    double loc[2] = {err, perr}, g[2];
    MPI_Allreduce(loc, g, 2, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0)
        std::printf("[E Mercury shell, np=%d] after %d steps: max |U - U_serial| = %.2e (scale %.2f), max |Phi - Phi_serial| = %.2e, "
                    "div B = %.2e (serial %.2e), mass_H = %.10e (serial %.10e)\n", np, nsteps, g[0], scale, g[1], d.div_b, dref.div_b, d.mass_H, dref.mass_H);
    CHECK(g[0] < 1e-9 * std::max(scale, 1.0));
    CHECK(g[1] < 1e-9);
    CHECK(d.div_b < 1e-13);
    CHECK(std::abs(d.mass_H - dref.mass_H) < 1e-9 * dref.mass_H);
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    test_freestream(MPI_COMM_WORLD, rank, np);
    test_wave(MPI_COMM_WORLD, rank, np);
    test_sod(MPI_COMM_WORLD, rank, np);
    test_alfven(MPI_COMM_WORLD, rank, np);
    test_shell(MPI_COMM_WORLD, rank, np);
    int fails = g_failures, all = 0;
    MPI_Allreduce(&fails, &all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", all ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return all ? 1 : 0;
}
