// mg_b_proto: two-level B multigrid prototype with frozen fluid (docs/Multigrid.md section 11).  Run on ONE rank.
//   mg_b_proto setup.txt mercury.txt coarse_dir=<dir with global.unsb> fluid_fields=<fields_final.unsb of the coarse mesh> mode=plain|mg
//              n_cycles=.. nu1=.. nuc=.. nu2=.. relax=.. report=..
// Fine mesh = uniform 2x2x2 octree refinement of the coarse mesh (spherical snap); the fluid state of the fine mesh is the coarse one prolonged piecewise constant
// and stays frozen.  B (induced flux) starts from zero on both levels.  The smoother is the full induction step of the solver (ideal + ambipolar EMF,
// implicit resistive shell / Hall) with ONE uniform step dt per level, so the fixed point of every level is d1 E = 0.
// Two-level FAS for the flux:  restriction = sum of the child face fluxes (exact discrete Stokes),  prolongation = transfer_B (d2 = 0 per fine cell).
#include <mpi.h>

#include <chrono>
#include <cmath>
#include <cstdio>

#include "driver/Solver.h"
#include "uns/core/Params.h"
#include "uns/io/MeshIO.h"
#include "uns/parallel/Distributor.h"
#include "uns/partition/Partitioner.h"
#include "uns/refine/Refine.h"

using namespace uns;
using namespace mercury;

namespace
{
LocalMesh serial_local(const Mesh &m)
{
    Geometry g = compute_geometry(m);
    CellGraph cg = build_cell_graph(m, g, CostModel());
    PartOptions po;
    po.method = PartMethod::Rcb;
    Distributor d(m, partition_cells(cg, 1, po), 1);
    return d.build(0);
}
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rc = 0;
    try
    {
        Params p;
        for (int i = 1; i < argc; ++i)
        {
            const std::string a = argv[i];
            if (a.find('=') == std::string::npos) p.load_file(a); else p.set_override(a);
        }
        Config cfg = Config::from_params(p);
        cfg.steady = true; cfg.mg_levels = 0; cfg.write_vtk = false; cfg.write_fields = false; cfg.write_profile = false; cfg.log_every = 0;
        const std::string cdir = p.get_string("coarse_dir", "geometry");
        const Mesh Cm = io::read_local_mesh(cdir + "/global.unsb").mesh;
        const Geometry cg = compute_geometry(Cm);
        RefineSpec spec;
        spec.max_level = 1; spec.solid_max_level = 1;
        RefineOptions opt;
        opt.spherical_snap = true;
        const RefineResult rr = build_refined_mesh(Cm, spec, std::vector<int>(Cm.cx.n_cells(), 1), opt);
        const Geometry fg = compute_geometry(rr.fine);
        std::printf("coarse %d cells / %d faces, fine %d cells / %d faces\n", Cm.cx.n_cells(), Cm.cx.n_faces(), rr.fine.cx.n_cells(), rr.fine.cx.n_faces());
        const LocalMesh lmC = serial_local(Cm), lmF = serial_local(rr.fine);
        Solver sC(MPI_COMM_WORLD, lmC, cfg), sF(MPI_COMM_WORLD, lmF, cfg);
        sC.initialize(); sF.initialize();
        // frozen fluid
        {
            std::map<std::string, io::FieldRecord> fl;
            for (auto &f : io::read_fields(p.get_string("fluid_fields", ""))) fl[f.name] = std::move(f);
            for (int k = 0; k < 2; ++k)
            {
                const std::vector<double> &U = fl[k == 0 ? "U_H" : "U_Na"].values;
                const CellComplex &cxc = lmC.mesh.cx, &cxf = lmF.mesh.cx;
                for (Id c = 0; c < cxc.n_cells(); ++c)
                    for (int m = 0; m < 5; ++m) sC.state().U[k][5 * static_cast<std::size_t>(c) + m] = U[5 * static_cast<std::size_t>(cxc.ent[3].gid[c]) + m];
                for (Id c = 0; c < cxf.n_cells(); ++c)
                    for (int m = 0; m < 5; ++m) sF.state().U[k][5 * static_cast<std::size_t>(c) + m] = U[5 * static_cast<std::size_t>(rr.parent_cell[cxf.ent[3].gid[c]]) + m];
            }
        }
        const CellComplex &cxC = lmC.mesh.cx, &cxF = lmF.mesh.cx;
        auto gather = [](const CellComplex &cx, const std::vector<double> &loc) { // local face vector -> global
            std::vector<double> g(cx.n_faces());
            for (Id f = 0; f < cx.n_faces(); ++f) g[cx.ent[2].gid[f]] = loc[f];
            return g;
        };
        auto scatter = [](const CellComplex &cx, const std::vector<double> &glob) {
            std::vector<double> l(cx.n_faces());
            for (Id f = 0; f < cx.n_faces(); ++f) l[f] = glob[cx.ent[2].gid[f]];
            return l;
        };
        std::vector<double> sgn(rr.fine.cx.n_faces(), 0.0);
        for (Id f = 0; f < rr.fine.cx.n_faces(); ++f)
            if (rr.face_coarse_face[f] >= 0) sgn[f] = dot(fg.face_area[f], cg.face_area[rr.face_coarse_face[f]]) > 0 ? 1.0 : -1.0;
        auto restrict_faces = [&](const std::vector<double> &fine_global) {
            std::vector<double> c(Cm.cx.n_faces(), 0.0);
            for (Id f = 0; f < rr.fine.cx.n_faces(); ++f)
                if (rr.face_coarse_face[f] >= 0) c[rr.face_coarse_face[f]] += sgn[f] * fine_global[f];
            return c;
        };
        const std::string mode = p.get_string("mode", "plain");
        const int ncyc = p.get_int("n_cycles", 100), nu1 = p.get_int("nu1", 4), nuc = p.get_int("nuc", 4), nu2 = p.get_int("nu2", 0), rep = p.get_int("report", 10);
        const double relax = p.get_double("relax", 1.0), dt_scale = p.get_double("dt_scale", 1.0);
        const double dtF = dt_scale * sF.b_dt(), dtC = dt_scale * sC.b_dt();
        std::printf("mode %s: dt_fine %.3e dt_coarse %.3e (ratio %.2f), nu1 %d nuc %d nu2 %d relax %.2f\n", mode.c_str(), dtF, dtC, dtC / dtF, nu1, nuc, nu2, relax);
        std::printf("%8s %10s %10s %12s\n", "fine_steps", "wall[s]", "work[F]", "res_B");
        const double t0 = now();
        long fine_steps = 0, coarse_steps = 0;
        double res = 0;
        for (int cyc = 1; cyc <= ncyc; ++cyc)
        {
            if (mode == "plain")
            {
                for (int s = 0; s < nu1 + nu2; ++s) { res = sF.b_smooth(dtF); ++fine_steps; }
            }
            else
            {
                for (int s = 0; s < nu1; ++s) { res = sF.b_smooth(dtF); ++fine_steps; }
                // FAS two-level correction
                std::vector<double> rf, rc;
                sF.b_rate(dtF, rf);
                const std::vector<double> phiF = gather(cxF, sF.state().Phi);
                const std::vector<double> phic0 = restrict_faces(phiF);
                const std::vector<double> rate_c = restrict_faces(gather(cxF, rf));
                sC.state().Phi = scatter(cxC, phic0);
                sC.b_rate(dtC, rc);
                std::vector<double> g(cxC.n_faces());
                const std::vector<double> rcg = gather(cxC, rc);
                std::vector<double> gg(Cm.cx.n_faces());
                for (std::size_t i = 0; i < gg.size(); ++i) gg[i] = rate_c[i] - rcg[i];
                g = scatter(cxC, gg);
                for (int s = 0; s < nuc; ++s) { sC.b_smooth(dtC, &g); ++coarse_steps; }
                std::vector<double> delta = gather(cxC, sC.state().Phi);
                for (std::size_t i = 0; i < delta.size(); ++i) delta[i] -= phic0[i];
                const std::vector<double> dfine = transfer_B(Cm, cg, delta, rr, fg);
                const std::vector<double> dl = scatter(cxF, dfine);
                for (Id f = 0; f < cxF.n_faces(); ++f) sF.state().Phi[f] += relax * dl[f];
                for (int s = 0; s < nu2; ++s) { res = sF.b_smooth(dtF); ++fine_steps; }
            }
            if (cyc % rep == 0 || cyc == ncyc)
                std::printf("%8ld %10.1f %10.1f %12.4e\n", fine_steps, now() - t0, fine_steps + (coarse_steps + 2.0 * cyc * (mode == "mg")) / 8.0, res);
        }
    }
    catch (const std::exception &e)
    {
        std::fprintf(stderr, "mg_b_proto: error: %s\n", e.what());
        rc = 1;
    }
    MPI_Finalize();
    return rc;
}
