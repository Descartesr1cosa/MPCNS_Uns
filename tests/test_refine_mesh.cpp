// Non-conforming refinement: topology, closure, volumes, hanging-node polyhedra, balance;
// B transfer (div-free per fine cell, exact for a uniform field); conservative U prolongation.
#include <cstdio>
#include <random>

#include "TestUtil.h"
#include "uns/complex/Validate.h"
#include "uns/generate/HexMesh.h"
#include "uns/refine/Refine.h"

using namespace uns;

static double total_volume(const Geometry &g)
{
    double v = 0;
    for (double x : g.cell_volume) v += x;
    return v;
}

static void check_mesh(const char *name, const Mesh &coarse, const RefineSpec &spec, const std::vector<int> &lev, bool planar)
{
    RefineResult rr = build_refined_mesh(coarse, spec, lev);
    const TopologyReport t = validate_topology(rr.fine.cx);
    const Geometry fg = compute_geometry(rr.fine), cg = compute_geometry(coarse);
    const GeometryReport gr = check_geometry(rr.fine, fg);
    long poly = 0, maxf = 0;
    for (Id c = 0; c < rr.fine.cx.n_cells(); ++c)
    {
        poly += rr.fine.cx.c2f.count(c) > 6;
        maxf = std::max<long>(maxf, rr.fine.cx.c2f.count(c));
    }
    std::printf("[%s] coarse=%d fine=%d (balance iters %d) polyhedra=%ld max_faces=%ld\n        %s\n        %s\n", name,
                coarse.cx.n_cells(), rr.fine.cx.n_cells(), rr.balance_iterations, poly, maxf, t.summary().c_str(),
                gr.summary().c_str());
    CHECK(t.ok);
    CHECK(t.euler == (planar ? 1 : 2));
    CHECK(gr.ok);
    CHECK(gr.misoriented_dual == 0);
    if (planar) CHECK(std::abs(total_volume(fg) - total_volume(cg)) < 1e-12 * total_volume(cg));
    // boundary tags survive
    long tagged_c = 0, tagged_f = 0;
    for (int tg : coarse.cx.face_tag) tagged_c += tg != 0;
    for (int tg : rr.fine.cx.face_tag) tagged_f += tg != 0;
    CHECK(tagged_f >= tagged_c);

    // ---- B transfer: B = d1(A) on the coarse mesh is divergence free ----
    std::mt19937 rng(3);
    std::vector<double> A(coarse.cx.n_edges()), Bc(coarse.cx.n_faces(), 0.0);
    for (double &a : A) a = std::uniform_real_distribution<double>(-1, 1)(rng);
    for (Id f = 0; f < coarse.cx.n_faces(); ++f)
        for (auto k = coarse.cx.f2e.begin(f); k < coarse.cx.f2e.end(f); ++k) Bc[f] += coarse.cx.f2e.sgn[k] * A[coarse.cx.f2e.idx[k]];
    TransferReport rep;
    const auto Bf = transfer_B(coarse, cg, Bc, rr, fg, &rep);
    std::printf("        B transfer: max rel |d2 B| = %.2e, max tile-sum defect = %.2e\n", rep.max_rel_div, rep.max_tile_defect);
    CHECK(rep.max_rel_div < 1e-13);
    CHECK(rep.max_tile_defect < 1e-13);
    (void)Bf;

    // ---- uniform field: coarse flux = b.A_f exactly (planar faces) -> fine flux must equal b.A_f ----
    if (planar)
    {
        const Vec3 b0(0.3, -1.1, 0.7);
        std::vector<double> Bu(coarse.cx.n_faces());
        for (Id f = 0; f < coarse.cx.n_faces(); ++f) Bu[f] = dot(b0, cg.face_area[f]);
        TransferReport r2;
        const auto Bfu = transfer_B(coarse, cg, Bu, rr, fg, &r2);
        double err = 0;
        for (Id f = 0; f < rr.fine.cx.n_faces(); ++f) err = std::max(err, std::abs(Bfu[f] - dot(b0, fg.face_area[f])));
        std::printf("        uniform field reproduced to %.2e\n", err);
        CHECK(err < 1e-12);
    }

    // ---- U prolongation: conservation and bounds ----
    std::vector<double> U(coarse.cx.n_cells());
    for (Id c = 0; c < coarse.cx.n_cells(); ++c) U[c] = 2.0 + std::sin(0.7 * cg.cell_center[c].x) + 0.3 * cg.cell_center[c].y;
    const auto Uf = prolong_cells(coarse, cg, U, 1, rr, fg, true);
    std::vector<double> mass(coarse.cx.n_cells(), 0.0), vol(coarse.cx.n_cells(), 0.0);
    for (Id c = 0; c < rr.fine.cx.n_cells(); ++c)
    {
        mass[rr.parent_cell[c]] += Uf[c] * fg.cell_volume[c];
        vol[rr.parent_cell[c]] += fg.cell_volume[c];
    }
    double cons = 0;
    for (Id c = 0; c < coarse.cx.n_cells(); ++c) cons = std::max(cons, std::abs(mass[c] / vol[c] - U[c]));
    std::printf("        U prolongation: max |mean(U_fine) - U_coarse| = %.2e\n", cons);
    CHECK(cons < 1e-12);
}

int main()
{
    BoxSpec bs;
    bs.nx = 4; bs.ny = 3; bs.nz = 3; bs.hi = {4, 3, 3};
    Mesh box = build_mesh(make_box(bs));
    {   // uniform level 2: 8^2 children per cell
        RefineSpec s; s.max_level = 3;
        std::vector<int> lev(box.cx.n_cells(), 2);
        RefineResult rr = build_refined_mesh(box, s, lev);
        CHECK(rr.fine.cx.n_cells() == box.cx.n_cells() * 64);
        check_mesh("box uniform L2", box, s, lev, true);
    }
    {   // one cell at level 3 in a level-0 sea: needs 2:1 balancing
        RefineSpec s; s.max_level = 3;
        std::vector<int> lev(box.cx.n_cells(), 0);
        lev[17] = 3;
        check_mesh("box single deep cell", box, s, lev, true);
    }
    {   // random levels
        RefineSpec s; s.max_level = 3;
        std::mt19937 rng(11);
        std::vector<int> lev(box.cx.n_cells());
        for (int &l : lev) l = rng() % 3;
        check_mesh("box random", box, s, lev, true);
    }
    {   // manual region cuts through coarse cells (leaf-level refinement inside a cell)
        RefineSpec s; s.max_level = 3;
        Region r; r.kind = Region::Kind::Sphere; r.a = {1.2, 1.3, 1.1, 0.8}; r.level = 3;
        s.regions.push_back(r);
        std::vector<int> lev(box.cx.n_cells(), 0);
        check_mesh("box sphere region", box, s, lev, true);
    }
    {   // curved shell with valence-3 corners, solid core, spherical snap
        CubedSphereSpec cs;
        cs.ncube = 3; cs.nr = 4; cs.n_solid = 1; cs.r_in = 1; cs.r_out = 5; cs.r_stretch = 1.2;
        Mesh shell = build_mesh(make_cubed_sphere_shell(cs));
        RefineSpec s; s.max_level = 2; s.solid_max_level = 0;
        Region r; r.kind = Region::Kind::Box; r.a = {-6, -6, 0.5, 6, 6, 6}; r.level = 2;
        s.regions.push_back(r);
        std::vector<int> lev(shell.cx.n_cells(), 0);
        RefineOptions o; o.spherical_snap = true;
        RefineResult rr = build_refined_mesh(shell, s, lev, o);
        const TopologyReport t = validate_topology(rr.fine.cx);
        const GeometryReport gr = check_geometry(rr.fine, compute_geometry(rr.fine));
        std::printf("[shell box region, snap] fine=%d  %s\n        %s\n", rr.fine.cx.n_cells(), t.summary().c_str(), gr.summary().c_str());
        CHECK(t.ok && t.euler == 2 && gr.ok);
        check_mesh("shell region", shell, s, lev, false);
    }
    TEST_END();
}
