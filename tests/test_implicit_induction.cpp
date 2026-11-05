// Semi-implicit induction (E as unknown): Hall whistler packet against the analytic dispersion, unconditional stability,
// exact div B for any solver accuracy, consistency with the explicit operator for small dt, and np-independence.
#include <mpi.h>

#include <cmath>
#include <complex>
#include <cstdio>

#include "TestUtil.h"
#include "diagnostics/DivB.h"
#include "induction/Faraday.h"
#include "induction/ImplicitInduction.h"
#include "uns/generate/HexMesh.h"
#include "uns/geometry/Reconstruct.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"

using namespace uns;
using namespace mercury;
using cplx = std::complex<double>;

namespace
{
const double kGamma = 0.5;   // alpha * B0 (hall_coef 0.5, n_e = 1, B0 = 1)
const double kK = 2 * 3.14159265358979 / 1.0, kSigma = 1.0;

LocalMesh local_of(const Mesh &m, int np, int rank)
{
    Geometry g = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, g, CostModel());
    PartOptions po;
    po.method = PartMethod::Rcb;
    Distributor d(m, partition_cells(cg, np, po), np);
    return d.build(rank);
}

// psi(x,t) of psi_t = -i gamma psi_xx for the initial Gaussian packet exp(-x^2/2s^2) exp(i k x), by quadrature in Fourier space
cplx psi_exact(double x, double t)
{
    cplx sum = 0;
    const int n = 4000;
    const double qmax = kK + 8.0 / kSigma, dq = 2 * qmax / n;
    for (int i = 0; i < n; ++i)
    {
        const double q = -qmax + (i + 0.5) * dq;
        const double amp = kSigma * std::exp(-0.5 * kSigma * kSigma * (q - kK) * (q - kK));
        sum += amp * std::exp(cplx(0, kGamma * q * q * t + q * x)) * dq;
    }
    return sum * std::sqrt(2 * 3.14159265358979) / (2 * 3.14159265358979);
}

struct Run
{
    std::vector<double> Phi;
    std::vector<double> by, bz; // cell centred perturbation at the end
    double divb = 0;
    int iters = 0, steps = 0;
    double max_res = 0;
    double norm0 = 0, norm1 = 0;
};

Run evolve(MPI_Comm comm, const LocalMesh &lm, double dt, double T, double rtol)
{
    const CellComplex &cx = lm.mesh.cx;
    const Geometry geo = compute_geometry(lm.mesh);
    Config c;
    c.hall = true;
    c.sc.hall_coef = 0.5;
    c.hall_r_min = -100; c.hall_r_max = -99; // taper = 1
    c.ind_rtol = rtol; c.ind_max_it = 400; c.ind_restart = 30;
    comm::Exchanger ex(comm, lm);
    ImplicitInduction ind(lm, geo, ex, c);
    std::vector<double> B0(3 * static_cast<std::size_t>(cx.n_cells()), 0.0), ne(cx.n_cells(), 1.0);
    for (Id i = 0; i < cx.n_cells(); ++i) B0[3 * i] = 1.0;
    ind.set_frozen_state(B0.data(), ne.data());

    Run r;
    r.Phi.assign(cx.n_faces(), 0.0);
    // psi0 = b_y + i b_z  (chi = x-dependence only, so div b = 0 face by face)
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const double x = geo.face_center[f].x;
        const cplx p = std::exp(-x * x / (2 * kSigma * kSigma)) * std::exp(cplx(0, kK * x));
        r.Phi[f] = dot(Vec3(0, p.real(), p.imag()), geo.face_area[f]);
    }
    std::vector<double> E(cx.n_edges(), 0.0), Eexp(cx.n_edges(), 0.0);
    std::vector<char> dir(cx.n_edges(), 0);
    const auto active = active_faces(cx);
    auto energy = [&] {
        double s = 0;
        for (Id f = 0; f < cx.ent[2].n_owned; ++f) s += r.Phi[f] * r.Phi[f] * geo.star2[f];
        ex.allreduce(&s, 1, MPI_SUM);
        return s;
    };
    r.norm0 = energy();
    const int nsteps = static_cast<int>(std::lround(T / dt));
    for (int s = 0; s < nsteps; ++s)
    {
        ex.forward(Dim::Face, r.Phi.data());
        double res = 0;
        r.iters += ind.solve(r.Phi.data(), Eexp.data(), dir, dt, E.data(), &res);
        r.max_res = std::max(r.max_res, res);
        faraday_update(cx, active, E.data(), r.Phi.data(), dt);
        ++r.steps;
    }
    ex.forward(Dim::Face, r.Phi.data());
    r.norm1 = energy();
    r.divb = div_b(cx, r.Phi.data()).max_rel;
    ex.allreduce(&r.divb, 1, MPI_MAX);
    const auto b = reconstruct_cell_vector(cx, geo, r.Phi.data(), cx.n_cells());
    r.by.resize(cx.n_cells()); r.bz.resize(cx.n_cells());
    for (Id i = 0; i < cx.n_cells(); ++i) { r.by[i] = b[i].y; r.bz[i] = b[i].z; }
    return r;
}

double l2_error(const LocalMesh &lm, const Run &r, double T, MPI_Comm comm)
{
    const CellComplex &cx = lm.mesh.cx;
    const Geometry geo = compute_geometry(lm.mesh);
    double e = 0, n = 0;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        const cplx ex = psi_exact(geo.cell_center[c].x, T);
        e += std::norm(cplx(r.by[c], r.bz[c]) - ex);
        n += std::norm(ex);
    }
    double v[2] = {e, n};
    MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_DOUBLE, MPI_SUM, comm);
    return std::sqrt(v[0] / v[1]);
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    BoxSpec bs;
    bs.nx = 200; bs.ny = 2; bs.nz = 2; bs.lo = {-8, 0, 0}; bs.hi = {8, 0.16, 0.16};
    const Mesh m = build_mesh(make_box(bs));
    const LocalMesh lm = local_of(m, np, rank);
    const double T = 0.1, h = 0.08;
    const double dt_explicit = h * h / (kGamma * 3.14159265358979 * 3.14159265358979); // whistler limit of an explicit scheme (order of magnitude)

    // (1) accuracy against the analytic whistler packet: space error (coarse vs fine mesh) and time error (dt, dt/2 on the fine mesh)
    BoxSpec bf = bs;
    bf.nx = 400;
    const Mesh mf = build_mesh(make_box(bf));
    const LocalMesh lmf = local_of(mf, np, rank);
    const Run a = evolve(MPI_COMM_WORLD, lm, 0.0025, T, 1e-10);
    const Run b = evolve(MPI_COMM_WORLD, lmf, 0.0025, T, 1e-10);
    const Run b2 = evolve(MPI_COMM_WORLD, lmf, 0.00125, T, 1e-10);
    const double ea = l2_error(lm, a, T, MPI_COMM_WORLD), eb = l2_error(lmf, b, T, MPI_COMM_WORLD), eb2 = l2_error(lmf, b2, T, MPI_COMM_WORLD);
    if (rank == 0)
        std::printf("[whistler packet, np=%d] L2 error, h=0.08: %.3e ; h=0.04: %.3e (space order %.2f) ; dt 0.0025 -> 0.00125 on h=0.04: %.3e -> %.3e ; "
                    "GMRES %d it/step, div B %.1e\n", np, ea, eb, std::log2(ea / eb), eb, eb2, a.iters / a.steps, a.divb);
    CHECK(ea < 0.3);
    CHECK(std::log2(ea / eb) > 1.2);   // converges with h (the observed order includes the first-order time error of the backward-Euler step)
    CHECK(eb2 < eb);                   // and the time error is smaller still on the finer step
    CHECK(a.divb < 1e-13 && b.divb < 1e-13 && b2.divb < 1e-13);
    CHECK(a.max_res < 1e-9);

    // (2) unconditional stability and exact div B with a loose solver: dt = 40 x the explicit limit, rtol 1e-3
    const Run c = evolve(MPI_COMM_WORLD, lm, 0.1, 0.5, 1e-3);
    if (rank == 0)
        std::printf("[large dt = %.0fx explicit limit, rtol 1e-3] energy %.4e -> %.4e (never grows), %d it/step, div B %.1e\n", 0.1 / dt_explicit, c.norm0, c.norm1,
                    c.iters / c.steps, c.divb);
    CHECK(c.norm1 <= c.norm0 * 1.0000001);
    CHECK(c.divb < 1e-13); // exact whatever the linear-solver accuracy

    // (3) np independence of the result of the accurate run
    {
        const LocalMesh ref = local_of(m, 1, 0);
        const Run s = evolve(MPI_COMM_SELF, ref, 0.0025, T, 1e-10);
        double err = 0, sc = 0;
        for (Id f = 0; f < lm.mesh.cx.ent[2].n_owned; ++f)
        {
            const GId g = lm.mesh.cx.ent[2].gid[f];
            err = std::max(err, std::abs(a.Phi[f] - s.Phi[g]));
            sc = std::max(sc, std::abs(s.Phi[g]));
        }
        MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::printf("[np independence] max |Phi_parallel - Phi_serial| = %.2e (scale %.3f)\n", err, sc);
        CHECK(err < 1e-8 * std::max(sc, 1e-30) + 1e-12);
    }
    int fails = g_failures, all = 0;
    MPI_Allreduce(&fails, &all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", all ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return all ? 1 : 0;
}
