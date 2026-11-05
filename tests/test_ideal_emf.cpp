// Ideal EMF on edges of arbitrary valence.
//  (1) regular hexahedra: bitwise-level agreement with a transcription of the legacy structured
//      Mercury scheme (AssembleOneDirectionEMF_ + AssembleEdgeEMF_FromFaceE_Ideal_) on random data
//  (2) uniform B: the transverse upwind part vanishes on every edge (valence 3, 5, hanging nodes, ...)
//  (3) uniform u, B on planar polyhedral meshes: E = -(u x B).dl exactly
//  (4) np-independence of E and of B after CT steps, d2 B at round-off (MPI, any number of ranks)
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <random>

#include "induction/Faraday.h"
#include "induction/IdealEMF.h"
#include "diagnostics/DivB.h"
#include "TestUtil.h"
#include "uns/generate/HexMesh.h"
#include "uns/geometry/Reconstruct.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"
#include "uns/refine/Refine.h"

using namespace uns;
using namespace mercury;

namespace
{
struct Setup
{
    LocalMesh lm;
    Geometry geo;
    std::vector<double> u, Bc, Phi;
};

LocalMesh serial_local(const Mesh &m)
{
    std::vector<int> part(m.cx.n_cells(), 0);
    Distributor d(m, part, 1);
    return d.build(0);
}

std::vector<double> d1_of(const CellComplex &cx, const std::vector<double> &A)
{
    std::vector<double> B(cx.n_faces(), 0.0);
    for (Id f = 0; f < cx.n_faces(); ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) B[f] += cx.f2e.sgn[k] * A[cx.f2e.idx[k]];
    return B;
}

// run the ideal EMF on the owned edges of lm and return E (owned edges valid)
std::vector<double> ideal_E(MPI_Comm comm, const LocalMesh &lm, const Geometry &geo, const std::vector<double> &u,
                            const std::vector<double> &Phi, double cuct, const std::vector<double> *Bc_override = nullptr)
{
    const CellComplex &cx = lm.mesh.cx;
    std::vector<double> Bc3(3 * static_cast<std::size_t>(cx.n_cells()));
    const auto Bv = reconstruct_cell_vector(cx, geo, Phi.data(), cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c) { Bc3[3 * c] = Bv[c].x; Bc3[3 * c + 1] = Bv[c].y; Bc3[3 * c + 2] = Bv[c].z; }
    if (Bc_override) Bc3 = *Bc_override;
    comm::Exchanger ex(comm, lm);
    EdgeEMF emf(lm, geo, ex);
    std::vector<double> E(cx.n_edges(), 0.0);
    emf.begin(E.data());
    IdealInputs in;
    in.u_cell = u.data(); in.B_cell_tot = Bc3.data(); in.B_cell_ind = Bc3.data(); in.phi_tot = Phi.data(); in.phi_ind = Phi.data();
    IdealParams p; p.cuct = cuct;
    add_ideal_emf(lm, geo, in, p, E.data());
    emf.finish(E.data());
    return E;
}

// n kites around one axis (edge valence n on the axis), extruded in z: planar faces, irregular valence
Mesh make_wheel(int n, int nz)
{
    HexMesh h;
    const double pi = 4 * std::atan(1.0);
    auto add_layer_nodes = [&](int k) {
        h.x.push_back({0, 0, 0.5 * k});
        for (int i = 0; i < n; ++i) h.x.push_back({std::cos(2 * pi * i / n), std::sin(2 * pi * i / n), 0.5 * k});
        for (int i = 0; i < n; ++i) h.x.push_back({1.7 * std::cos(2 * pi * (i + 0.5) / n), 1.7 * std::sin(2 * pi * (i + 0.5) / n), 0.5 * k});
    };
    for (int k = 0; k <= nz; ++k) add_layer_nodes(k);
    const int per = 1 + 2 * n;
    for (int k = 0; k < nz; ++k)
        for (int i = 0; i < n; ++i)
        {
            auto id = [&](int kk, int which) { return static_cast<Id>(kk * per + which); }; // 0: O, 1+i: P_i, 1+n+i: Q_i
            const int p0 = 1 + i, p1 = 1 + (i + 1) % n, q = 1 + n + i;
            std::array<Id, 8> v = {id(k, 0), id(k, p0), id(k, p1), id(k, q), id(k + 1, 0), id(k + 1, p0), id(k + 1, p1), id(k + 1, q)};
            h.hex.push_back(v);
            h.material.push_back(kFluid);
        }
    return build_mesh(h);
}

bool interior_edge(const CellComplex &cx, const Geometry &geo, Id e)
{
    return cx.e2f.count(e) == static_cast<int>(geo.ec_ptr[e + 1] - geo.ec_ptr[e]);
}

// ---------------------------------------------------------------- (1) legacy transcription
void test_legacy_regression()
{
    const int n = 6;
    const double h = 0.5;
    BoxSpec bs;
    bs.nx = bs.ny = bs.nz = n; bs.hi = {n * h, n * h, n * h};
    Mesh m = build_mesh(make_box(bs));
    LocalMesh lm = serial_local(m);
    Geometry geo = compute_geometry(lm.mesh);
    const CellComplex &cx = lm.mesh.cx;
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(-1, 1);
    std::vector<double> A(cx.n_edges());
    for (double &a : A) a = U(rng);
    std::vector<double> Phi = d1_of(cx, A);
    std::vector<double> u(3 * cx.n_cells());
    for (double &x : u) x = U(rng);

    // legacy arrays: ucon, flux per axis at (i,j,k)
    auto idx = [&](int i, int j, int k) { return (k * n + j) * n + i; };
    std::vector<std::array<double, 3>> uc(n * n * n);
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        const Vec3 x = geo.cell_center[c];
        const int i = int(std::floor(x.x / h)), j = int(std::floor(x.y / h)), k = int(std::floor(x.z / h));
        for (int a = 0; a < 3; ++a) uc[idx(i, j, k)][a] = u[3 * c + a] / h;
    }
    std::vector<double> Fx((n + 1) * n * n), Fy(n * (n + 1) * n), Fz(n * n * (n + 1));
    auto fx = [&](int i, int j, int k) -> double & { return Fx[(k * n + j) * (n + 1) + i]; };
    auto fy = [&](int i, int j, int k) -> double & { return Fy[(k * (n + 1) + j) * n + i]; };
    auto fz = [&](int i, int j, int k) -> double & { return Fz[(k * n + j) * n + i]; };
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        const Vec3 c = geo.face_center[f], a = geo.face_area[f];
        const int ax = std::abs(a.x) > std::abs(a.y) && std::abs(a.x) > std::abs(a.z) ? 0 : (std::abs(a.y) > std::abs(a.z) ? 1 : 2);
        const double sgn = (ax == 0 ? a.x : ax == 1 ? a.y : a.z) > 0 ? 1.0 : -1.0;
        const int i = int(std::floor(c.x / h + 1e-9)), j = int(std::floor(c.y / h + 1e-9)), k = int(std::floor(c.z / h + 1e-9));
        (ax == 0 ? fx(i, j, k) : ax == 1 ? fy(i, j, k) : fz(i, j, k)) = sgn * Phi[f];
    }
    auto Fl = [&](int axis, int i, int j, int k) { return axis == 0 ? fx(i, j, k) : axis == 1 ? fy(i, j, k) : fz(i, j, k); };
    auto ucon = [&](int a, int i, int j, int k) { return uc[idx(i, j, k)][a]; };
    auto bind = [&](int a, int i, int j, int k) {
        return 0.5 * (Fl(a, i, j, k) + Fl(a, i + (a == 0), j + (a == 1), k + (a == 2)));
    };
    auto G = [&](int beta, int fax, int iL, int jL, int kL, int iR, int jR, int kR, double bface) {
        const double GL = ucon(fax, iL, jL, kL) * bind(beta, iL, jL, kL) - ucon(beta, iL, jL, kL) * bface;
        const double GR = ucon(fax, iR, jR, kR) * bind(beta, iR, jR, kR) - ucon(beta, iR, jR, kR) * bface;
        const double rad = std::max(std::abs(ucon(fax, iL, jL, kL)), std::abs(ucon(fax, iR, jR, kR)));
        return 0.5 * (GL + GR) - 0.5 * rad * (bind(beta, iR, jR, kR) - bind(beta, iL, jL, kL));
    };
    // face candidates e[3] at face (axis,i,j,k)
    auto Eface = [&](int axis, int i, int j, int k, int comp) -> double {
        const int iL = i - (axis == 0), jL = j - (axis == 1), kL = k - (axis == 2);
        const double bf = Fl(axis, i, j, k);
        double e[3] = {0, 0, 0};
        if (axis == 0) { e[1] = G(2, 0, iL, jL, kL, i, j, k, bf); e[2] = -G(1, 0, iL, jL, kL, i, j, k, bf); }
        if (axis == 1) { e[0] = -G(2, 1, iL, jL, kL, i, j, k, bf); e[2] = G(0, 1, iL, jL, kL, i, j, k, bf); }
        if (axis == 2) { e[0] = G(1, 2, iL, jL, kL, i, j, k, bf); e[1] = -G(0, 2, iL, jL, kL, i, j, k, bf); }
        return e[comp];
    };
    auto amax = [&](int a, int i0, int j0, int k0, int i1, int j1, int k1, int i2, int j2, int k2, int i3, int j3, int k3) {
        return std::max(std::max(std::abs(ucon(a, i0, j0, k0)), std::abs(ucon(a, i1, j1, k1))),
                        std::max(std::abs(ucon(a, i2, j2, k2)), std::abs(ucon(a, i3, j3, k3))));
    };
    const double C = 0.5;
    auto legacy = [&](int axis, int i, int j, int k) -> double {
        if (axis == 0)
        {
            const double Ec = 0.25 * (Eface(1, i, j, k, 0) + Eface(1, i, j, k - 1, 0) + Eface(2, i, j, k, 0) + Eface(2, i, j - 1, k, 0));
            const double ae = amax(1, i, j, k, i, j - 1, k, i, j, k - 1, i, j - 1, k - 1);
            const double az = amax(2, i, j, k, i, j - 1, k, i, j, k - 1, i, j - 1, k - 1);
            return Ec + C * (0.5 * ae * (fz(i, j, k) - fz(i, j - 1, k)) - 0.5 * az * (fy(i, j, k) - fy(i, j, k - 1)));
        }
        if (axis == 1)
        {
            const double Ec = 0.25 * (Eface(0, i, j, k, 1) + Eface(0, i, j, k - 1, 1) + Eface(2, i, j, k, 1) + Eface(2, i - 1, j, k, 1));
            const double az = amax(2, i, j, k, i - 1, j, k, i, j, k - 1, i - 1, j, k - 1);
            const double ax = amax(0, i, j, k, i - 1, j, k, i, j, k - 1, i - 1, j, k - 1);
            return Ec + C * (0.5 * az * (fx(i, j, k) - fx(i, j, k - 1)) - 0.5 * ax * (fz(i, j, k) - fz(i - 1, j, k)));
        }
        const double Ec = 0.25 * (Eface(0, i, j, k, 2) + Eface(0, i, j - 1, k, 2) + Eface(1, i, j, k, 2) + Eface(1, i - 1, j, k, 2));
        const double ax = amax(0, i, j, k, i - 1, j, k, i, j - 1, k, i - 1, j - 1, k);
        const double ae = amax(1, i, j, k, i - 1, j, k, i, j - 1, k, i - 1, j - 1, k);
        return Ec + C * (0.5 * ax * (fy(i, j, k) - fy(i - 1, j, k)) - 0.5 * ae * (fx(i, j, k) - fx(i, j - 1, k)));
    };

    const std::vector<double> E = ideal_E(MPI_COMM_SELF, lm, geo, u, Phi, C);
    double maxerr = 0, maxE = 0;
    long compared = 0;
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        const Vec3 t = geo.edge_vec[e], mid = geo.edge_mid[e];
        const int ax = std::abs(t.x) > 0.1 ? 0 : (std::abs(t.y) > 0.1 ? 1 : 2);
        const int i = int(std::floor(mid.x / h + 1e-9)), j = int(std::floor(mid.y / h + 1e-9)), k = int(std::floor(mid.z / h + 1e-9));
        // all four surrounding cells / faces must exist
        const int c1 = (ax == 0 ? j : i), c2 = (ax == 2 ? j : k);
        if (ax == 0 && !(j >= 1 && j <= n - 1 && k >= 1 && k <= n - 1)) continue;
        if (ax == 1 && !(i >= 1 && i <= n - 1 && k >= 1 && k <= n - 1)) continue;
        if (ax == 2 && !(i >= 1 && i <= n - 1 && j >= 1 && j <= n - 1)) continue;
        (void)c1; (void)c2;
        const double ref = legacy(ax, i, j, k);
        maxerr = std::max(maxerr, std::abs(E[e] - ref));
        maxE = std::max(maxE, std::abs(ref));
        ++compared;
    }
    std::printf("[legacy regression] %ld interior edges, max |E_new - E_legacy| = %.3e (max |E| = %.3f)\n", compared, maxerr, maxE);
    CHECK(compared > 100);
    CHECK(maxerr < 1e-13 * std::max(maxE, 1.0));
}

// ---------------------------------------------------------------- (2),(3) uniform fields
double test_uniform(const char *name, const Mesh &m, bool planar)
{
    LocalMesh lm = serial_local(m);
    Geometry geo = compute_geometry(lm.mesh);
    const CellComplex &cx = lm.mesh.cx;
    const Vec3 b0(0.4, -0.9, 0.6), u0(0.7, 0.2, -0.5);
    std::vector<double> Phi(cx.n_faces());
    for (Id f = 0; f < cx.n_faces(); ++f) Phi[f] = dot(b0, geo.face_area[f]);
    std::mt19937 rng(9);
    std::uniform_real_distribution<double> U(-1, 1);
    std::vector<double> ur(3 * cx.n_cells()), u1(3 * cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c)
        for (int a = 0; a < 3; ++a) { ur[3 * c + a] = U(rng); u1[3 * c + a] = (a == 0 ? u0.x : a == 1 ? u0.y : u0.z); }
    std::vector<double> Bc(3 * cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c) { Bc[3 * c] = b0.x; Bc[3 * c + 1] = b0.y; Bc[3 * c + 2] = b0.z; }

    // (2) random u, uniform B: UCT part must vanish -> E(C=0.5) == E(C=0)
    const auto Ea = ideal_E(MPI_COMM_SELF, lm, geo, ur, Phi, 0.5, &Bc);
    const auto Eb = ideal_E(MPI_COMM_SELF, lm, geo, ur, Phi, 0.0, &Bc);
    double d = 0, s = 0, dbnd = 0;
    int vmax = 0, vmin = 99;
    long nint = 0;
    std::map<int, long> valence;
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        if (!interior_edge(cx, geo, e)) { dbnd = std::max(dbnd, std::abs(Ea[e] - Eb[e])); continue; }
        ++nint;
        ++valence[cx.e2f.count(e)];
        vmax = std::max(vmax, cx.e2f.count(e));
        vmin = std::min(vmin, cx.e2f.count(e));
        d = std::max(d, std::abs(Ea[e] - Eb[e]));
        s = std::max(s, std::abs(Ea[e]));
    }
    std::printf("[%s] %ld interior edges, valence %d..%d {", name, nint, vmin, vmax);
    for (auto &kv : valence) std::printf(" %d:%ld", kv.first, kv.second);
    std::printf(" }; uniform B, random u: |UCT part| max = %.2e (E scale %.2f); boundary edges: %.2e\n", d, s, dbnd);
    // planar faces: exact; warped faces: the face normals have O(warp) components along the edge, the residual is a
    // consistent first-order effect (a * |B| * warp) and must stay small
    CHECK(d < (planar ? 1e-12 : 8e-2) * std::max(s, 1.0));

    // (3) uniform u, B: E = -(u x B).dl
    const auto Ec = ideal_E(MPI_COMM_SELF, lm, geo, u1, Phi, 0.5, &Bc);
    double err = 0, sc = 0, mean_rel = 0;
    long cnt = 0;
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        if (!interior_edge(cx, geo, e)) continue;
        const double ex = -dot(cross(u0, b0), geo.edge_vec[e]);
        err = std::max(err, std::abs(Ec[e] - ex));
        sc = std::max(sc, std::abs(ex));
        mean_rel += std::abs(Ec[e] - ex) / (geo.edge_len[e] * norm(cross(u0, b0)));
        ++cnt;
    }
    mean_rel /= std::max(cnt, 1L);
    std::printf("[%s] uniform u,B: max |E + (u x B).dl| = %.2e (scale %.3f), mean relative error %.2e%s\n", name, err, sc, mean_rel, planar ? "" : "  (warped faces: first-order in the face warp)");
    if (planar) CHECK(err < 1e-12 * std::max(sc, 1.0));
    else CHECK(err < 0.08 * sc);
    return mean_rel;
}

// ---------------------------------------------------------------- (4) partition independence + CT
std::vector<double> ct_run(MPI_Comm comm, const LocalMesh &lm, const Geometry &geo, int steps, double dt, double *divmax)
{
    const CellComplex &cx = lm.mesh.cx;
    // smooth analytic data evaluated per entity (identical on every rank)
    std::vector<double> A(cx.n_edges());
    for (Id e = 0; e < cx.n_edges(); ++e)
    {
        const Vec3 x = geo.edge_mid[e];
        const Vec3 a(std::sin(0.9 * x.y + 0.3), std::cos(0.7 * x.z), std::sin(0.5 * x.x + 0.2 * x.y));
        A[e] = dot(a, geo.edge_vec[e]);
    }
    std::vector<double> Phi = d1_of(cx, A), u(3 * cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        const Vec3 x = geo.cell_center[c];
        u[3 * c] = 0.5 + 0.3 * std::sin(x.y); u[3 * c + 1] = 0.2 * std::cos(x.x + x.z); u[3 * c + 2] = -0.1 + 0.2 * std::sin(x.x);
    }
    comm::Exchanger ex(comm, lm);
    EdgeEMF emf(lm, geo, ex);
    const auto active = active_faces(cx);
    std::vector<double> E(cx.n_edges()), Bc3(3 * static_cast<std::size_t>(cx.n_cells()));
    *divmax = 0;
    for (int s = 0; s < steps; ++s)
    {
        ex.forward(Dim::Face, Phi.data());
        const auto Bv = reconstruct_cell_vector(cx, geo, Phi.data(), cx.n_cells());
        for (Id c = 0; c < cx.n_cells(); ++c) { Bc3[3 * c] = Bv[c].x; Bc3[3 * c + 1] = Bv[c].y; Bc3[3 * c + 2] = Bv[c].z; }
        emf.begin(E.data());
        IdealInputs in;
        in.u_cell = u.data(); in.B_cell_tot = Bc3.data(); in.B_cell_ind = Bc3.data(); in.phi_tot = Phi.data(); in.phi_ind = Phi.data();
        add_ideal_emf(lm, geo, in, IdealParams(), E.data());
        emf.finish(E.data());
        faraday_update(cx, active, E.data(), Phi.data(), dt);
        *divmax = std::max(*divmax, div_b(cx, Phi.data()).max_rel);
    }
    ex.forward(Dim::Face, Phi.data());
    return Phi;
}

void test_parallel(const char *name, const Mesh &m, MPI_Comm comm)
{
    int rank, np;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);
    // serial reference on every rank
    LocalMesh ref = serial_local(m);
    Geometry rg = compute_geometry(ref.mesh);
    double dref;
    const auto Pref = ct_run(MPI_COMM_SELF, ref, rg, 20, 0.02, &dref);

    Geometry mg = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, mg, CostModel());
    PartOptions po; po.method = PartMethod::Rcb;
    Distributor dist(m, partition_cells(cg, np, po), np);
    LocalMesh lm = dist.build(rank);
    Geometry geo = compute_geometry(lm.mesh);
    double dloc;
    const auto P = ct_run(comm, lm, geo, 20, 0.02, &dloc);
    double err = 0, sc = 0;
    for (Id f = 0; f < lm.mesh.cx.n_faces(); ++f)
    {
        const GId g = lm.mesh.cx.ent[2].gid[f];
        err = std::max(err, std::abs(P[f] - Pref[g]));
        sc = std::max(sc, std::abs(Pref[g]));
    }
    double gerr = 0, gdiv = 0;
    MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&dloc, &gdiv, 1, MPI_DOUBLE, MPI_MAX, comm);
    if (rank == 0)
        std::printf("[%s] np=%d: max |B_parallel - B_serial| after 20 steps = %.2e (scale %.2f), max rel d2B = %.2e (serial %.2e)\n",
                    name, np, gerr, sc, gdiv, dref);
    CHECK(gerr < 1e-12 * std::max(sc, 1.0));
    CHECK(gdiv < 1e-13);
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    BoxSpec bs;
    bs.nx = 5; bs.ny = 4; bs.nz = 4; bs.hi = {5, 4, 4};
    Mesh box = build_mesh(make_box(bs));
    // hanging-node polyhedra with planar faces
    RefineSpec rs; rs.max_level = 2;
    Region reg; reg.kind = Region::Kind::Sphere; reg.a = {2.2, 1.8, 2.1, 1.3}; reg.level = 2;
    rs.regions.push_back(reg);
    Mesh hang = build_refined_mesh(box, rs, std::vector<int>(box.cx.n_cells(), 0)).fine;
    CubedSphereSpec cs;
    cs.ncube = 4; cs.nr = 4; cs.n_solid = 0; // all fluid: edges of solid-only cells are skipped by the ideal EMF cs.r_in = 1; cs.r_out = 5; cs.r_stretch = 1.2;
    Mesh shell = build_mesh(make_cubed_sphere_shell(cs));

    if (rank == 0)
    {
        test_legacy_regression();
        test_uniform("box", box, true);
        test_uniform("hanging", hang, true);
        double c[3];
        int k = 0;
        for (int nc : {4, 8, 16})
        {
            CubedSphereSpec q = cs;
            q.ncube = nc; q.nr = nc; q.r_stretch = 1.0;
            c[k++] = test_uniform(("shell " + std::to_string(nc)).c_str(), build_mesh(make_cubed_sphere_shell(q)), false);
        }
        std::printf("[shell] mean relative error of uniform u x B: %.2e -> %.2e -> %.2e on halving h\n", c[0], c[1], c[2]);
        CHECK(c[1] < 0.65 * c[0] && c[2] < 0.65 * c[1]); // consistent: first order in the face warp
        for (int nw : {3, 5, 6, 8}) test_uniform(("wheel" + std::to_string(nw)).c_str(), make_wheel(nw, 3), true);
    }
    test_parallel("hanging", hang, MPI_COMM_WORLD);
    test_parallel("shell", shell, MPI_COMM_WORLD);
    test_parallel("wheel5", make_wheel(5, 6), MPI_COMM_WORLD);
    test_parallel("wheel6", make_wheel(6, 6), MPI_COMM_WORLD);
    int fails = g_failures, all = 0;
    MPI_Allreduce(&fails, &all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::printf("%s\n", all ? "TEST FAILED" : "TEST PASSED");
    MPI_Finalize();
    return all ? 1 : 0;
}
