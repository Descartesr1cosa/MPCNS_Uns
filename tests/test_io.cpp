// Round trip of the binary container (sub-complex + plans + fields) and a VTU file
// that tools/check_vtu.py validates (polyhedron volumes vs. cell volumes).
#include <cstdio>
#include <cstdlib>

#include "TestUtil.h"
#include "uns/complex/Validate.h"
#include "uns/generate/HexMesh.h"
#include "uns/io/MeshIO.h"
#include "uns/io/VtkWriter.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"

using namespace uns;

int main(int argc, char **argv)
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    CubedSphereSpec cs;
    cs.ncube = 4; cs.nr = 4; cs.n_solid = 1; cs.r_in = 0.8; cs.r_out = 5.0; cs.r_stretch = 1.2;
    Mesh m = build_mesh(make_cubed_sphere_shell(cs));
    Geometry geo = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, geo, CostModel());
    const int np = 3;
    auto part = partition_cells(cg, np, PartOptions());
    Distributor dist(m, part, np);
    for (int r = 0; r < np; ++r)
    {
        LocalMesh a = dist.build(r);
        const std::string f = dir + "/test_rank" + std::to_string(r) + ".unsb";
        io::write_local_mesh(f, a);
        LocalMesh b = io::read_local_mesh(f);
        CHECK(topology_hash(a.mesh.cx) == topology_hash(b.mesh.cx));
        CHECK(a.mesh.x.size() == b.mesh.x.size());
        CHECK(a.mesh.cx.ent[1].gid == b.mesh.cx.ent[1].gid);
        CHECK(a.plan.reduce_send.idx == b.plan.reduce_send.idx);
        CHECK(a.plan.recv[3].peers == b.plan.recv[3].peers);
        CHECK(a.mesh.tag_names == b.mesh.tag_names);
        CHECK(validate_topology(b.mesh.cx).ok);
        // piece for ParaView, with rank id and cell volume as cell data
        Geometry g = compute_geometry(b.mesh);
        io::VtkField fr{"rank", 1, std::vector<double>(b.n_owned(Dim::Cell), r)};
        io::VtkField fv{"volume", 1, std::vector<double>(g.cell_volume.begin(), g.cell_volume.begin() + b.n_owned(Dim::Cell))};
        io::VtkField fm{"material", 1, {}};
        for (Id c = 0; c < b.n_owned(Dim::Cell); ++c) fm.values.push_back(b.mesh.cx.cell_material[c]);
        io::write_vtu(dir + "/test_rank" + std::to_string(r) + ".vtu", b.mesh, b.n_owned(Dim::Cell), {fr, fv, fm});
    }
    io::write_pvtu(dir + "/test.pvtu", {"test_rank0.vtu", "test_rank1.vtu", "test_rank2.vtu"},
                   {{"rank", 1, {}}, {"volume", 1, {}}, {"material", 1, {}}});
    io::write_pvd(dir + "/test.pvd", {{0.0, "test.pvtu"}});

    std::vector<io::FieldRecord> fs = {{"B", Dim::Face, 1, std::vector<double>(m.cx.n_faces(), 0.5)},
                                       {"U", Dim::Cell, 5, std::vector<double>(m.cx.n_cells() * 5, 1.5)}};
    io::write_fields(dir + "/test_fields.unsb", topology_hash(m.cx), fs);
    std::uint64_t h = 0;
    auto back = io::read_fields(dir + "/test_fields.unsb", &h);
    CHECK(h == topology_hash(m.cx));
    CHECK(back.size() == 2 && back[1].ncomp == 5 && back[1].values.size() == fs[1].values.size());
    TEST_END();
}
