#include "Pipeline.h"

#include <sys/stat.h>

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

using namespace uns;

namespace pre
{
namespace
{
void make_dir(const std::string &d) { ::mkdir(d.c_str(), 0755); }

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

int run_generate(const Params &p)
{
    const std::string out = p.get_string("output_dir", "geometry");
    const int np = p.get_int("proc_num", 4);
    make_dir(out);
    Mesh m = make_mesh(p);
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
    text += "\nrank   owned(N,E,F,C)                        ghost(N,E,F,C)                      halo send (cells)  edge-reduce send\n";
    std::vector<std::string> pieces;
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
        text += b;
        if (p.get_bool("write_vtk", true))
        {
            write_piece_vtu(out, r, lm, part);
            std::snprintf(b, sizeof b, "part_%05d.vtu", r);
            pieces.push_back(b);
        }
    }
    if (!pieces.empty())
    {
        io::write_pvtu(out + "/partition.pvtu", pieces, {{"rank", 1, {}}, {"volume", 1, {}}, {"material", 1, {}}});
    }
    // global complex (used by the coarse-solution / refinement stages)
    LocalMesh gm;
    gm.mesh = m;
    io::write_local_mesh(out + "/global.unsb", gm);
    std::ofstream(out + "/partition_report.txt") << text;
    std::printf("%s", text.c_str());
    std::printf("wrote %d sub-complexes to %s\n", np, out.c_str());
    return 0;
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
