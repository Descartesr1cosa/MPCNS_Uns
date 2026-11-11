#include "app/App.h"

#ifdef _WIN32
#include <process.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "diagnostics/DivB.h"
#include "induction/EdgeEMF.h"
#include "induction/Faraday.h"
#include "induction/HallEMF.h"
#include "state/FieldSet.h"
#include "uns/complex/Validate.h"
#include "uns/io/MeshIO.h"

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
    mhd::FieldSet fs(cx);
    auto &B = fs.add("B", Dim::Face);
    auto &E = fs.add("E", Dim::Edge);
    auto &A = fs.add("A", Dim::Edge);
    for (Id e = 0; e < cx.n_edges(); ++e) A.data[e] = synthetic_emf(7, cx.ent[1].gid[e]);
    for (Id f = 0; f < cx.n_faces(); ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) B.data[f] += cx.f2e.sgn[k] * A.data[cx.f2e.idx[k]];
    const auto active = mhd::active_faces(cx);
    mhd::EdgeEMF emf(lm, geo, ex);
    mhd::HallConfig hall;
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
        mhd::add_hall_emf(hall, emf, E.ptr()); // no-op until the Hall term exists
        emf.finish(E.ptr());                       // owner reduction + broadcast
        mhd::faraday_update(cx, active, E.ptr(), B.ptr(), 1.0);
        const auto d = mhd::div_b(cx, B.ptr());
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

namespace
{
bool under_mpi()
{
    for (const char *v : {"OMPI_COMM_WORLD_SIZE", "PMIX_RANK", "PMI_SIZE", "PMI_RANK", "MPI_LOCALRANKID", "MV2_COMM_WORLD_SIZE", "SLURM_PROCID", "UNS_RELAUNCHED"})
        if (std::getenv(v)) return true;
    return false;
}

std::string find_in_path(const std::string &prog)
{
    namespace fs = std::filesystem;
    if (prog.find('/') != std::string::npos || prog.find('\\') != std::string::npos) return fs::exists(prog) ? prog : "";
    const char *path = std::getenv("PATH");
    std::istringstream is(path ? path : "");
#ifdef _WIN32
    const char sep = ';';
    const char *exts[] = {"", ".exe", ".bat", ".cmd"};
#else
    const char sep = ':';
    const char *exts[] = {""};
#endif
    for (std::string d; std::getline(is, d, sep);)
        for (const char *e : exts)
            if (!d.empty() && fs::exists(fs::path(d) / (prog + e))) return (fs::path(d) / (prog + e)).string();
    return "";
}

// `./MercuryUns` started without mpirun: read config.yml, work out the number of ranks (nproc, or "auto" = number of
// part_*.unsb in geometry_dir) and re-execute ourselves under the MPI launcher.  Returns only when no relaunch is needed.
void maybe_relaunch(int argc, char **argv)
{
    if (under_mpi()) return;
    namespace fs = std::filesystem;
    const fs::path cwd = fs::current_path();
    Params p;
    try
    {
        if (!p.load_args(std::vector<std::string>(argv + 1, argv + argc))) return;
    }
    catch (const std::exception &) { fs::current_path(cwd); return; } // the real parse below reports the error
    int np = 1;
    const std::string n = p.get_string("nproc", "auto");
    if (n == "auto")
    {
        np = 0;
        std::error_code ec;
        for (const auto &e : fs::directory_iterator(p.get_string("geometry_dir", "geometry"), ec))
        {
            const std::string f = e.path().filename().string();
            if (f.rfind("part_", 0) == 0 && f.size() > 5 && f.substr(f.size() - 5) == ".unsb") ++np;
        }
        if (np == 0) np = 1;
    }
    else np = std::max(1, std::atoi(n.c_str()));
    fs::current_path(cwd); // the child parses the same arguments again from the original directory
    if (np == 1 && p.get_string("mpi_launcher", "").empty()) return; // single rank: MPI singleton start
    #ifdef _WIN32
    std::string launcher = p.get_string("mpi_launcher", "mpiexec"); // MS-MPI
#else
    std::string launcher = p.get_string("mpi_launcher", "mpirun");
#endif
    std::string exe = find_in_path(launcher);
#ifdef UNS_MPIEXEC
    if (exe.empty() && (launcher == "mpirun" || launcher == "mpiexec")) exe = UNS_MPIEXEC; // the one this binary was built with
#endif
    if (exe.empty())
    {
        std::fprintf(stderr, "%s: cannot find the MPI launcher '%s'; load your MPI module (PATH) or set mpi_launcher in config.yml\n", argv[0], launcher.c_str());
        std::exit(1);
    }
    std::vector<std::string> a{exe};
    const std::string base = fs::path(exe).filename().string();
    if (p.get_string("mpi_np_flag", "").empty()) a.push_back("-np"); else a.push_back(p.get_string("mpi_np_flag", ""));
    a.push_back(std::to_string(np));
#ifndef _WIN32
    if (geteuid() == 0 && (base == "mpirun" || base == "mpiexec") && p.get_string("mpi_args", "").find("allow-run-as-root") == std::string::npos &&
        fs::exists(fs::path(exe).parent_path() / "ompi_info")) // Open MPI refuses to run as root without this
        a.push_back("--allow-run-as-root");
#endif
    {
        std::istringstream is(p.get_string("mpi_args", ""));
        for (std::string t; is >> t;) a.push_back(t);
    }
    {
        std::error_code ec;
        fs::path self = fs::read_symlink("/proc/self/exe", ec); // Linux
#ifdef _WIN32
        char buf[MAX_PATH * 4];
        const DWORD n = GetModuleFileNameA(nullptr, buf, sizeof buf);
        if (n > 0 && n < sizeof buf) self = std::string(buf, n);
#endif
        a.push_back(self.empty() ? fs::absolute(argv[0]).string() : self.string());
    }
    for (int i = 1; i < argc; ++i) a.push_back(argv[i]);
#ifdef _WIN32
    _putenv_s("UNS_RELAUNCHED", "1");
#else
    setenv("UNS_RELAUNCHED", "1", 1);
#endif
    std::printf("[launch] %s\n", [&] { std::string s; for (auto &t : a) s += t + " "; return s; }().c_str());
    std::fflush(stdout);
    std::vector<char *> cargv;
    for (auto &t : a) cargv.push_back(t.data());
    cargv.push_back(nullptr);
#ifdef _WIN32
    // no exec() on Windows: run the launcher, wait, pass its exit code on
    const intptr_t rc = _spawnv(_P_WAIT, exe.c_str(), cargv.data());
    if (rc == -1) std::perror("spawn");
    std::exit(rc == -1 ? 1 : static_cast<int>(rc));
#else
    execv(exe.c_str(), cargv.data());
    std::perror("execv");
    std::exit(1);
#endif
}
} // namespace

namespace mhd
{
int run_app(int argc, char **argv, const char *name, const RunFn &run)
{
    maybe_relaunch(argc, argv);
    MPI_Init(&argc, &argv);
    int rank = 0, np = 1, rc = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    try
    {
        Params p;
        // config.yml (default) and/or setup files (later ones win), then name=value overrides; relative paths are relative to the yml's folder
        if (!p.load_args(std::vector<std::string>(argv + 1, argv + argc)))
            throw std::runtime_error("no config.yml in the current directory and no setup file given");
        char fname[64];
        std::snprintf(fname, sizeof fname, "/part_%05d.unsb", rank);
        LocalMesh lm = io::read_local_mesh(p.get_string("geometry_dir", "geometry") + fname);
        if (lm.nranks != np) throw std::runtime_error("geometry was partitioned for " + std::to_string(lm.nranks) + " ranks, running on " + std::to_string(np));
        const TopologyReport t = validate_topology(lm.mesh.cx);
        if (!t.ok) throw std::runtime_error("invalid local complex: " + t.summary());
        const std::string task = p.get_string("task", "selftest");
        if (task == "selftest") rc = selftest(MPI_COMM_WORLD, lm, p);
        else if (task == "run") run(MPI_COMM_WORLD, lm, p);
        else throw std::runtime_error("task '" + task + "' is not implemented yet (selftest | run)");
    }
    catch (const std::exception &e)
    {
        std::fprintf(stderr, "[rank %d] %s: error: %s\n", rank, name, e.what());
        rc = 1;
    }
    int grc = rc;
    MPI_Allreduce(&rc, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return grc;
}
} // namespace mhd
