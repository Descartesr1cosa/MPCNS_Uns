#include "Pipeline.h"

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "uns/complex/Validate.h"
#include "uns/generate/HexMesh.h"
#include "uns/geometry/Geometry.h"
#include "uns/io/MeshIO.h"
#include "uns/io/VtkWriter.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"
#include "uns/refine/Refine.h"
#include <map>

using namespace uns;

namespace pre
{
namespace
{
void make_dir(const std::string &d) // mkdir -p
{
    for (std::size_t i = 1; i <= d.size(); ++i)
        if (i == d.size() || d[i] == '/') ::mkdir(d.substr(0, i).c_str(), 0755);
}

std::string rank_file(const std::string &dir, int r)
{
    char b[64];
    std::snprintf(b, sizeof b, "/part_%05d.unsb", r);
    return dir + b;
}

Mesh make_mesh(const Params &p)
{
    const std::string src = p.get_string("mesh_source", "box");
    if (src == "box")
    {
        BoxSpec s;
        s.nx = p.get_int("box_nx", 16); s.ny = p.get_int("box_ny", 16); s.nz = p.get_int("box_nz", 16);
        const auto lo = p.get_doubles("box_lo"), hi = p.get_doubles("box_hi");
        if (lo.size() == 3) s.lo = {lo[0], lo[1], lo[2]};
        if (hi.size() == 3) s.hi = {hi[0], hi[1], hi[2]};
        return build_mesh(make_box(s));
    }
    if (src == "cubed_sphere")
    {
        CubedSphereSpec s;
        s.ncube = p.get_int("cs_ncube", 8); s.nr = p.get_int("cs_nr", 8); s.n_solid = p.get_int("cs_nsolid", 0);
        s.r_in = p.get_double("cs_r_in", 1.0); s.r_out = p.get_double("cs_r_out", 10.0);
        s.r_stretch = p.get_double("cs_r_stretch", 1.0);
        // planet + ellipsoidal domain (Mercury): cs_r_in..cs_r_surf solid (uniform), r_surf..ellipsoid fluid (geometric)
        s.r_surf = p.get_double("cs_r_surf", 0.0);
        s.ellipsoid = p.get_bool("cs_ellipsoid", false);
        const auto lo = p.get_doubles("cs_ell_lo"), hi = p.get_doubles("cs_ell_hi");
        if (lo.size() == 3) s.ell_lo = {lo[0], lo[1], lo[2]};
        if (hi.size() == 3) s.ell_hi = {hi[0], hi[1], hi[2]};
        return build_mesh(make_cubed_sphere_shell(s));
    }
    throw std::runtime_error("unknown mesh_source '" + src + "' (box|cubed_sphere; plot3d/gmsh readers: see PLAN.md)");
}

CostModel cost_from(const Params &p)
{
    CostModel c;
    c.fluid = p.get_double("cost_fluid", c.fluid);
    c.solid = p.get_double("cost_solid", c.solid);
    c.per_extra_face = p.get_double("cost_extra_face", c.per_extra_face);
    c.comm_fluid = p.get_double("comm_fluid", c.comm_fluid);
    c.comm_solid = p.get_double("comm_solid", c.comm_solid);
    return c;
}

void write_piece_vtu(const std::string &dir, int r, const LocalMesh &lm, const std::vector<int> &part)
{
    (void)part;
    const Geometry g = compute_geometry(lm.mesh);
    const Id n = lm.n_owned(Dim::Cell);
    io::VtkField fr{"rank", 1, std::vector<double>(n, r)};
    io::VtkField fv{"volume", 1, std::vector<double>(g.cell_volume.begin(), g.cell_volume.begin() + n)};
    io::VtkField fm{"material", 1, {}};
    for (Id c = 0; c < n; ++c) fm.values.push_back(lm.mesh.cx.cell_material[c]);
    char b[64];
    std::snprintf(b, sizeof b, "/part_%05d.vtu", r);
    io::write_vtu(dir + b, lm.mesh, n, {fr, fv, fm});
}
} // namespace

namespace
{
// Validate, partition and write a (global) mesh as per-rank sub-complexes.
// `cell_fields` are global cell arrays (indexed by global cell id) added to the ParaView pieces.
int emit_case(const Mesh &m, const Params &p, const std::string &out, const std::vector<io::VtkField> &cell_fields)
{
    const int np = p.get_int("proc_num", 4);
    make_dir(out);
    const TopologyReport t = validate_topology(m.cx);
    std::printf("[topology] %s\n", t.summary().c_str());
    if (!t.ok) return 1;
    const Geometry geo = compute_geometry(m);
    const GeometryReport gr = check_geometry(m, geo);
    std::printf("[geometry] %s\n", gr.summary().c_str());
    if (!gr.ok) return 1;

    const CellGraph cg = build_cell_graph(m, geo, cost_from(p));
    PartOptions po;
    po.method = parse_method(p.get_string("partition_method", "hilbert"));
    po.refine_passes = p.get_int("partition_refine_passes", po.refine_passes);
    po.tolerance = p.get_double("partition_tolerance", po.tolerance);
    const std::vector<int> part = partition_cells(cg, np, po);
    const PartitionReport rep = evaluate_partition(cg, part, np);
    std::string text = rep.format(to_string(po.method));

    Distributor dist(m, part, np);
    text += "\nrank   owned(N,E,F,C)                        ghost(N,E,F,C)                      halo send (cells)  edge-reduce send   fluid   solid\n";
    std::vector<std::string> pieces;
    std::vector<io::VtkField> proto = {{"rank", 1, {}}, {"volume", 1, {}}, {"material", 1, {}}};
    for (const auto &f : cell_fields) proto.push_back({f.name, f.ncomp, {}});
    for (int r = 0; r < np; ++r)
    {
        LocalMesh lm = dist.build(r);
        io::write_local_mesh(rank_file(out, r), lm);
        char b[256];
        const auto &cx = lm.mesh.cx;
        std::snprintf(b, sizeof b, "%4d   %7d %7d %7d %7d     %7d %7d %7d %7d        %9zu        %9zu\n", r, cx.ent[0].n_owned,
                      cx.ent[1].n_owned, cx.ent[2].n_owned, cx.ent[3].n_owned, cx.count(Dim::Node) - cx.ent[0].n_owned,
                      cx.count(Dim::Edge) - cx.ent[1].n_owned, cx.count(Dim::Face) - cx.ent[2].n_owned,
                      cx.count(Dim::Cell) - cx.ent[3].n_owned, lm.plan.send[3].total(), lm.plan.reduce_send.total());
        {
            long nf = 0, ns = 0;
            for (Id c = 0; c < cx.ent[3].n_owned; ++c) (cx.cell_material[c] == kSolid ? ns : nf)++;
            char tail[64];
            std::snprintf(tail, sizeof tail, "  %7ld %7ld\n", nf, ns);
            std::string row(b);
            row.pop_back(); // the newline
            text += row + tail;
        }
        if (p.get_bool("write_vtk", true))
        {
            const Geometry g = compute_geometry(lm.mesh);
            const Id n = lm.n_owned(Dim::Cell);
            std::vector<io::VtkField> fs = {{"rank", 1, std::vector<double>(n, r)},
                                            {"volume", 1, std::vector<double>(g.cell_volume.begin(), g.cell_volume.begin() + n)},
                                            {"material", 1, {}}};
            for (Id c = 0; c < n; ++c) fs[2].values.push_back(cx.cell_material[c]);
            for (const auto &gf : cell_fields)
            {
                io::VtkField f{gf.name, gf.ncomp, {}};
                for (Id c = 0; c < n; ++c)
                    for (int k = 0; k < gf.ncomp; ++k) f.values.push_back(gf.values[static_cast<std::size_t>(cx.ent[3].gid[c]) * gf.ncomp + k]);
                fs.push_back(std::move(f));
            }
            std::snprintf(b, sizeof b, "/part_%05d.vtu", r);
            io::write_vtu(out + b, lm.mesh, n, fs);
            std::snprintf(b, sizeof b, "part_%05d.vtu", r);
            pieces.push_back(b);
        }
    }
    if (!pieces.empty()) io::write_pvtu(out + "/partition.pvtu", pieces, proto);
    LocalMesh gm; // global complex: input of the coarse-solution / refinement stages
    gm.mesh = m;
    io::write_local_mesh(out + "/global.unsb", gm);
    std::ofstream(out + "/partition_report.txt") << text;
    std::printf("%s", text.c_str());
    std::printf("wrote %d sub-complexes to %s\n", np, out.c_str());
    return 0;
}
} // namespace

int run_generate(const Params &p)
{
    return emit_case(make_mesh(p), p, p.get_string("output_dir", "geometry"), {});
}

// Synthetic coarse solution (dipole-like divergence-free B = d1(A), density bump in the tail, current sheet)
// used to exercise the refinement pipeline before the solver can produce a real one.
int run_demo_fields(const Params &p)
{
    const std::string dir = p.get_string("coarse_dir", "geometry");
    const Mesh m = io::read_local_mesh(dir + "/global.unsb").mesh;
    const Geometry g = compute_geometry(m);
    const CellComplex &cx = m.cx;
    const Vec3 mom(0, 0, -1.0);
    auto A_of = [&](const Vec3 &x) { // vector potential of a dipole + uniform field
        const double r = std::max(norm(x), 1e-12);
        return cross(mom, x) * (1.0 / (r * r * r)) + 0.5 * cross(Vec3(0, 0, 0.05), x);
    };
    std::vector<double> A(cx.n_edges()), B(cx.n_faces(), 0.0);
    for (Id e = 0; e < cx.n_edges(); ++e) A[e] = dot(A_of(g.edge_mid[e]), g.edge_vec[e]);
    for (Id f = 0; f < cx.n_faces(); ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) B[f] += cx.f2e.sgn[k] * A[cx.f2e.idx[k]];
    const Id nc = cx.n_cells();
    std::vector<double> rho(nc), pr(nc), Bc(3 * nc), Jc(3 * nc);
    for (Id c = 0; c < nc; ++c)
    {
        const Vec3 x = g.cell_center[c];
        rho[c] = 1.0 + 4.0 * std::exp(-((x.x + 3) * (x.x + 3) / 4.0 + x.y * x.y / 2.0 + x.z * x.z / 0.5));
        pr[c] = rho[c];
        const double sheet = std::exp(-(x.z * x.z) / 0.2) * (x.x < -1.5 ? 1.0 : 0.0);
        Bc[3 * c] = 0; Bc[3 * c + 1] = 0; Bc[3 * c + 2] = 1.0;
        Jc[3 * c] = 0; Jc[3 * c + 1] = 3.0 * sheet; Jc[3 * c + 2] = 0;
    }
    io::write_fields(p.get_string("coarse_fields", dir + "/coarse_fields.unsb"), topology_hash(cx),
                     {{"B", Dim::Face, 1, B}, {"rho", Dim::Cell, 1, rho}, {"p", Dim::Cell, 1, pr},
                      {"B_cell", Dim::Cell, 3, Bc}, {"J_cell", Dim::Cell, 3, Jc}});
    std::printf("wrote synthetic coarse fields (B, rho, p, B_cell, J_cell)\n");
    return 0;
}

// Coarse complex (+ optional coarse solution) -> refined non-conforming complex (+ transferred B and cell fields).
int run_refine(const Params &p)
{
    const std::string cdir = p.get_string("coarse_dir", "geometry");
    const std::string out = p.get_string("output_dir", "fine/geometry");
    const Mesh coarse = io::read_local_mesh(cdir + "/global.unsb").mesh;
    const Geometry cg = compute_geometry(coarse);
    const RefineSpec spec = RefineSpec::from_params(p);

    std::map<std::string, io::FieldRecord> fields;
    CoarseSolution sol;
    const std::string ffile = p.get_string("coarse_fields", "");
    if (!ffile.empty())
    {
        std::uint64_t h = 0;
        for (auto &f : io::read_fields(ffile, &h)) fields[f.name] = std::move(f);
        if (h != 0 && h != topology_hash(coarse.cx)) throw  // h == 0: written by the solver, which cannot hash the global complex
             std::runtime_error("coarse fields do not belong to the coarse complex (topology hash differs)");
        if (fields.count("rho")) sol.rho = fields["rho"].values;
        if (fields.count("p")) sol.p = fields["p"].values;
        auto vec3 = [&](const char *n, std::vector<Vec3> &v) {
            if (!fields.count(n)) return;
            const auto &d = fields[n].values;
            for (std::size_t i = 0; i + 2 < d.size(); i += 3) v.push_back({d[i], d[i + 1], d[i + 2]});
        };
        vec3("B_cell", sol.B_cell);
        vec3("J_cell", sol.J_cell);
    }
    RefineOptions opt;
    opt.spherical_snap = p.get_bool("refine_spherical_snap", false);
    // Cell budget: the flow-indicator thresholds are relaxed (x1.5) until the fine mesh fits.
    // Manual regions are never relaxed; if they alone exceed the budget a warning is printed.
    RefineSpec cur = spec;
    RefineResult rr;
    Id previous = -1;
    for (int attempt = 0;; ++attempt)
    {
        const std::vector<int> levels = compute_target_levels(coarse, cg, cur, ffile.empty() ? nullptr : &sol);
        rr = build_refined_mesh(coarse, cur, levels, opt);
        if (spec.budget_cells <= 0 || rr.fine.cx.n_cells() <= spec.budget_cells) break;
        if (cur.indicators.empty() || attempt >= 12 || rr.fine.cx.n_cells() == previous)
        {
            std::fprintf(stderr, "[refine] warning: %d fine cells exceed refine_budget_cells=%ld and cannot be reduced further\n",
                         rr.fine.cx.n_cells(), spec.budget_cells);
            break;
        }
        previous = rr.fine.cx.n_cells();
        for (auto &r : cur.indicators) r.threshold *= 1.5;
        std::printf("[refine] %d cells > budget %ld: indicator thresholds x1.5 (attempt %d)\n", rr.fine.cx.n_cells(), spec.budget_cells, attempt + 1);
    }
    const Geometry fg = compute_geometry(rr.fine);

    std::map<int, long> per_level;
    long poly = 0;
    for (Id c = 0; c < rr.fine.cx.n_cells(); ++c)
    {
        ++per_level[rr.cell_level[c]];
        poly += rr.fine.cx.c2f.count(c) > 6;
    }
    std::printf("[refine] coarse cells %d -> fine cells %d (balance iterations %d), polyhedra %ld\n", coarse.cx.n_cells(),
                rr.fine.cx.n_cells(), rr.balance_iterations, poly);
    for (const auto &kv : per_level) std::printf("         level %d: %ld cells\n", kv.first, kv.second);

    std::vector<io::FieldRecord> out_fields;
    std::vector<io::VtkField> vtk;
    vtk.push_back({"level", 1, std::vector<double>(rr.cell_level.begin(), rr.cell_level.end())});
    vtk.push_back({"coarse_cell", 1, std::vector<double>(rr.parent_cell.begin(), rr.parent_cell.end())});
    if (fields.count("B"))
    {
        TransferReport tr;
        out_fields.push_back({"B", Dim::Face, 1, transfer_B(coarse, cg, fields["B"].values, rr, fg, &tr)});
        std::printf("[transfer B] max relative |d2 B| over fine cells = %.2e, tile-flux defect = %.2e\n", tr.max_rel_div, tr.max_tile_defect);
        if (tr.max_rel_div > 1e-9) { std::fprintf(stderr, "B transfer is not divergence free\n"); return 1; }
    }
    if (fields.count("_time")) out_fields.push_back(fields["_time"]); // restart time of the coarse run
    for (const auto &kv : fields)
    {
        if (kv.second.dim != Dim::Cell) continue;
        io::FieldRecord f{kv.first, Dim::Cell, kv.second.ncomp, {}};
        if (p.get_bool("refine_prolong_linear", false))
            f.values = prolong_cells(coarse, cg, kv.second.values, kv.second.ncomp, rr, fg, true); // limited linear, conservative
        else // piecewise constant: positive and conservative by construction (conservative variables may leave the physical set under linear prolongation)
        {
            f.values.resize(static_cast<std::size_t>(rr.fine.cx.n_cells()) * kv.second.ncomp);
            for (Id c = 0; c < rr.fine.cx.n_cells(); ++c)
                for (int m = 0; m < kv.second.ncomp; ++m)
                    f.values[static_cast<std::size_t>(c) * kv.second.ncomp + m] = kv.second.values[static_cast<std::size_t>(rr.parent_cell[c]) * kv.second.ncomp + m];
        }
        vtk.push_back({f.name, f.ncomp, f.values});
        out_fields.push_back(std::move(f));
    }
    const int rc = emit_case(rr.fine, p, out, vtk);
    if (rc == 0 && !out_fields.empty())
    {
        io::write_fields(out + "/fine_fields.unsb", topology_hash(rr.fine.cx), out_fields);
        std::printf("wrote transferred fields to %s/fine_fields.unsb\n", out.c_str());
    }
    return rc;
}

int run_inspect(const Params &p)
{
    const std::string dir = p.get_string("output_dir", "geometry");
    const int np = p.get_int("proc_num", 1);
    int bad = 0;
    for (int r = 0; r < np; ++r)
    {
        LocalMesh lm = io::read_local_mesh(rank_file(dir, r));
        const TopologyReport t = validate_topology(lm.mesh.cx);
        const Geometry g = compute_geometry(lm.mesh);
        const GeometryReport gr = check_geometry(lm.mesh, g);
        std::printf("rank %d: %s\n        %s\n", r, t.summary().c_str(), gr.summary().c_str());
        bad += !t.ok;
    }
    return bad ? 1 : 0;
}

int run_export(const Params &p)
{
    const std::string dir = p.get_string("output_dir", "geometry");
    const int np = p.get_int("proc_num", 1);
    std::vector<std::string> pieces;
    for (int r = 0; r < np; ++r)
    {
        write_piece_vtu(dir, r, io::read_local_mesh(rank_file(dir, r)), {});
        char b[64];
        std::snprintf(b, sizeof b, "part_%05d.vtu", r);
        pieces.push_back(b);
    }
    io::write_pvtu(dir + "/partition.pvtu", pieces, {{"rank", 1, {}}, {"volume", 1, {}}, {"material", 1, {}}});
    return 0;
}
} // namespace pre
