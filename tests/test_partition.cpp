// All partitioners: every rank non-empty, balance within tolerance, cut statistics sane.
#include <algorithm>
#include <cstdio>

#include "TestUtil.h"
#include "uns/generate/HexMesh.h"
#include "uns/partition/Partitioner.h"

using namespace uns;

int main()
{
    CubedSphereSpec cs;
    cs.ncube = 10; cs.nr = 12; cs.n_solid = 3; cs.r_in = 0.8; cs.r_out = 8.0; cs.r_stretch = 1.1;
    Mesh m = build_mesh(make_cubed_sphere_shell(cs));
    Geometry geo = compute_geometry(m);
    CellGraph g = build_cell_graph(m, geo, CostModel());
    for (PartMethod pm : {PartMethod::Morton, PartMethod::Hilbert, PartMethod::Rcb})
        for (int np : {1, 2, 7, 16})
        {
            PartOptions o;
            o.method = pm;
            auto part = partition_cells(g, np, o);
            auto rep = evaluate_partition(g, part, np);
            std::printf("%-8s np=%2d  imbalance=%.4f  cut_faces=%ld  max_nb=%d\n", to_string(pm).c_str(), np, rep.imbalance,
                        rep.total_cut_faces, *std::max_element(rep.neighbours.begin(), rep.neighbours.end()));
            for (int p = 0; p < np; ++p) CHECK(rep.cells[p] > 0);
            CHECK(rep.imbalance < 1.10);
            if (np == 1) CHECK(rep.total_cut_faces == 0);
        }
    // refinement must not make the cut worse
    PartOptions a, b;
    a.method = b.method = PartMethod::Morton;
    a.refine_passes = 0;
    auto ra = evaluate_partition(g, partition_cells(g, 16, a), 16);
    auto rb = evaluate_partition(g, partition_cells(g, 16, b), 16);
    std::printf("morton np=16 cut: raw=%.1f refined=%.1f\n", ra.total_cut_weight, rb.total_cut_weight);
    CHECK(rb.total_cut_weight <= ra.total_cut_weight);
    std::printf("%s", rb.format("morton+refine").c_str());
    TEST_END();
}
