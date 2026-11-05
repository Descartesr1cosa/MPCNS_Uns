// MercuryUns <setup-file | -> [name=value ...]     (run with mpirun, one rank per part file)
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>

#include "diagnostics/DivB.h"
#include "hall/HallEMF.h"
#include "induction/EdgeEMF.h"
#include "induction/Faraday.h"
#include "state/FieldSet.h"
#include "uns/core/Params.h"
#include "uns/io/MeshIO.h"
#include "uns/complex/Validate.h"

using namespace uns;

namespace
{
// Deterministic pseudo-random line-integral contribution of cell c to edge e.
double synthetic_emf(GId gc, GId ge) { return std::sin(0.37 * static_cast<double>(gc) + 1.91 * static_cast<double>(ge)); }

// Structural self-test of the induction path: B0 = d1(A) is divergence free; random EMF
// contributions are assembled through the owner and applied with B -= dt d1(E).
int selftest(MPI_Comm comm, const LocalMesh &lm, const Params &p)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    const CellComplex &cx = lm.mesh.cx;
    const Geometry geo = compute_geometry(lm.mesh);
    comm::Exchanger ex(comm, lm);
    mercury::FieldSet fs(cx);
    auto &B = fs.add("B", Dim::Face);
    auto &E = fs.add("E", Dim::Edge);
    auto &A = fs.add("A", Dim::Edge);
    for (Id e = 0; e < cx.n_edges(); ++e) A.data[e] = synthetic_emf(7, cx.ent[1].gid[e]);
    for (Id f = 0; f < cx.n_faces(); ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) B.data[f] += cx.f2e.sgn[k] * A.data[cx.f2e.idx[k]];
    const auto active = mercury::active_faces(cx);
    mercury::EdgeEMF emf(lm, geo, ex);
    mercury::HallConfig hall;
    hall.enabled = p.get_bool("hall_enabled", false);
    hall.coef = p.get_double("hall_coef", 0.0);

    const int nsteps = p.get_int("selftest_steps", 200);
    double worst = 0;
    for (int s = 0; s < nsteps; ++s)
    {
        emf.begin(E.ptr());
        emf.for_each_owned_cell_edge([&](Id c, Id e, double w) {
            // stand-in for the ideal EMF: weights from the dual-face shares, random-like values
            emf.add(E.ptr(), e, 1e-2 * w * synthetic_emf(cx.ent[3].gid[c] + s, cx.ent[1].gid[e]));
        });
        mercury::add_hall_emf(hall, emf, E.ptr()); // no-op until the Hall term exists
        emf.finish(E.ptr());                       // owner reduction + broadcast
        mercury::faraday_update(cx, active, E.ptr(), B.ptr(), 1.0);
        const auto d = mercury::div_b(cx, B.ptr());
        worst = std::max(worst, d.max_rel);
    }
    double gworst = 0;
    MPI_Allreduce(&worst, &gworst, 1, MPI_DOUBLE, MPI_MAX, comm);
    // partition-independent checksum of B (owned faces, weighted by global id)
    double chk = 0, gchk = 0;
    for (Id f = 0; f < cx.ent[2].n_owned; ++f) chk += B.data[f] * static_cast<double>(cx.ent[2].gid[f] % 97 + 1);
    MPI_Allreduce(&chk, &gchk, 1, MPI_DOUBLE, MPI_SUM, comm);
    if (rank == 0)
        std::printf("[selftest] %d steps, max relative |d2 B| over all cells and steps = %.3e   checksum(B) = %.12e\n", nsteps,
                    gworst, gchk);
    return gworst < 1e-12 ? 0 : 1;
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, np = 1, rc = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    try
    {
        Params p;
        if (argc > 1 && std::string(argv[1]) != "-") p.load_file(argv[1]);
        for (int i = 2; i < argc; ++i) p.set_override(argv[i]);
        char name[64];
        std::snprintf(name, sizeof name, "/part_%05d.unsb", rank);
        LocalMesh lm = io::read_local_mesh(p.get_string("geometry_dir", "geometry") + name);
        if (lm.nranks != np) throw std::runtime_error("geometry was partitioned for " + std::to_string(lm.nranks) + " ranks, running on " + std::to_string(np));
        const TopologyReport t = validate_topology(lm.mesh.cx);
        if (!t.ok) throw std::runtime_error("invalid local complex: " + t.summary());
        const std::string task = p.get_string("task", "selftest");
        if (task == "selftest") rc = selftest(MPI_COMM_WORLD, lm, p);
        else throw std::runtime_error("task '" + task + "' is not implemented yet (see PLAN.md, milestone M4+)");
    }
    catch (const std::exception &e)
    {
        std::fprintf(stderr, "[rank %d] MercuryUns: error: %s\n", rank, e.what());
        rc = 1;
    }
    int grc = rc;
    MPI_Allreduce(&rc, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return grc;
}
