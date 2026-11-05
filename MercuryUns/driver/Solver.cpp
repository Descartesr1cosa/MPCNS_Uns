#include "driver/Solver.h"

#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "diagnostics/DivB.h"
#include "induction/Current.h"
#include "induction/Faraday.h"
#include "physics/Gyration.h"
#include "initial/Initial.h"
#include "io/GlobalFields.h"
#include "uns/geometry/Reconstruct.h"
#include "uns/io/MeshIO.h"
#include "uns/io/VtkWriter.h"

namespace mercury
{
using namespace uns;

namespace
{
void make_dir(const std::string &d)
{
    for (std::size_t i = 1; i <= d.size(); ++i)
        if (i == d.size() || d[i] == '/') ::mkdir(d.substr(0, i).c_str(), 0755);
}
} // namespace

Solver::Solver(MPI_Comm comm, const LocalMesh &lm, const Config &cfg)
    : comm_(comm), lm_(lm), cfg_(cfg), geo_(compute_geometry(lm.mesh)), ex_(comm, lm), emf_(lm, geo_, ex_), amb_(lm, geo_), res_(lm, geo_, ex_, cfg), ind_(lm, geo_, ex_, cfg),
      bc_(BoundaryInfo::build(lm, geo_, cfg)), stencil_(lm, geo_, bc_), active_(active_faces(lm.mesh.cx))
{
    MPI_Comm_rank(comm, &rank_);
    fs_[0] = {cfg.fsH.rho, cfg.u_inf.x, cfg.u_inf.y, cfg.u_inf.z, cfg.fsH.p};
    fs_[1] = {cfg.fsNa.rho, cfg.u_Na().x, cfg.u_Na().y, cfg.u_Na().z, cfg.fsNa.p};
}

void Solver::initialize()
{
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    set_initial_state(lm_, geo_, cfg_, s_.U[0], s_.U[1], s_.Phi, &time_);
    set_background(lm_, geo_, cfg_, PhiAdd_, BaddCell_);
    s1_ = s2_ = s_;
    for (int k = 0; k < 2; ++k)
    {
        W_[k].assign(5 * nc, 0.0);
        grad_[k].assign(15 * nc, 0.0);
        r0_.dU[k].assign(5 * nc, 0.0);
        r1_.dU[k].assign(5 * nc, 0.0);
    }
    uplus_.assign(3 * nc, 0.0);
    Bind_.assign(3 * nc, 0.0);
    Btot_.assign(3 * nc, 0.0);
    J_.assign(3 * nc, 0.0);
    PhiTot_.assign(cx.n_faces(), 0.0);
    dpe_.assign(3 * nc, 0.0);
    photo_.assign(nc, 0.0);
    vst_.assign(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c) // static chemistry fields at the cell centres
    {
        const Vec3 x = geo_.cell_center[c];
        vst_[c] = in_shadow(x) ? cfg_.chem.sk2 : cfg_.chem.sk1;
        if (cx.cell_material[c] == kFluid) photo_[c] = photo_rate(x, cfg_.sc.L_ref_m, cfg_.chem);
    }
    r0_.E.assign(cx.n_edges(), 0.0);
    r1_.E.assign(cx.n_edges(), 0.0);
    nstep_ = 0;
    initialized_ = true;
}

void Solver::prepare(State &s)
{
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    ex_.forward(Dim::Cell, s.U[0].data(), 5);
    ex_.forward(Dim::Cell, s.U[1].data(), 5);
    ex_.forward(Dim::Face, s.Phi.data(), 1);
    for (int k = 0; k < 2; ++k)
        for (std::size_t c = 0; c < nc; ++c)
        {
            const Prim w = to_prim(&s.U[k][5 * c], cfg_.sc.gamma, cfg_.rho_floor, cfg_.p_floor);
            double* o = &W_[k][5 * c];
            o[0] = w.rho; o[1] = w.u; o[2] = w.v; o[3] = w.w; o[4] = w.p;
        }
    const auto bi = reconstruct_cell_vector(cx, geo_, s.Phi.data(), cx.n_cells());
    for (std::size_t c = 0; c < nc; ++c)
    {
        const double b[3] = {bi[c].x, bi[c].y, bi[c].z};
        for (int a = 0; a < 3; ++a) { Bind_[3 * c + a] = b[a]; Btot_[3 * c + a] = b[a] + BaddCell_[3 * c + a]; }
        // charge-weighted ion velocity; solid cells carry no flow
        double up[3] = {0, 0, 0};
        if (cx.cell_material[c] == kFluid)
        {
            const NumInfo ni = cfg_.sc.numinfo(W_[0][5 * c], W_[1][5 * c]);
            for (int a = 0; a < 3; ++a) up[a] = ni.chiH * W_[0][5 * c + 1 + a] + ni.chiNa * W_[1][5 * c + 1 + a];
        }
        for (int a = 0; a < 3; ++a) uplus_[3 * c + a] = up[a];
    }
    for (Id f = 0; f < cx.n_faces(); ++f) PhiTot_[f] = s.Phi[f] + PhiAdd_[f];
    cell_current(lm_, geo_, Bind_.data(), J_.data());
    if (cfg_.ambipolar) electron_pressure_gradient(lm_, geo_, W_[0].data(), W_[1].data(), dpe_.data());
}

double Solver::local_dt() const
{
    return local_time_step(lm_, geo_, cfg_, W_[0].data(), W_[1].data(), Btot_.data());
}

void Solver::rhs(State &s, Rhs &r)
{
    const CellComplex &cx = lm_.mesh.cx;
    for (int k = 0; k < 2; ++k)
    {
        if (cfg_.second_order)
        {
            stencil_.limited_gradients(lm_, geo_, bc_, W_[k].data(), k, fs_[k], true, grad_[k].data());
            ex_.forward(Dim::Cell, grad_[k].data(), 15);
        }
        else
            std::fill(grad_[k].begin(), grad_[k].end(), 0.0); // piecewise constant (first order)
        fluid_flux_rhs(lm_, geo_, bc_, cfg_, k, fs_[k], W_[k].data(), grad_[k].data(), r.dU[k].data());
    }
    SourceFields sf;
    sf.photo = photo_.data(); sf.vst = vst_.data(); sf.dpe = dpe_.data();
    for (int k = 0; k < 2; ++k)
        add_sources(lm_, cfg_, k, W_[0].data(), W_[1].data(), uplus_.data(), Btot_.data(), J_.data(), sf, dt_cur_, r.dU[k].data());

    // EMF of the owned edges: ideal part now, the Hall term plugs in here later (hall/HallEMF.h)
    emf_.begin(r.E.data());
    IdealInputs in;
    in.u_cell = uplus_.data();
    in.B_cell_tot = Btot_.data();
    in.B_cell_ind = Bind_.data();
    in.phi_tot = PhiTot_.data();
    in.phi_ind = s.Phi.data();
    add_ideal_emf(lm_, geo_, in, cfg_.ideal, r.E.data());
    if (cfg_.ambipolar) amb_.add(cfg_, W_[0].data(), W_[1].data(), r.E.data());
    const Vec3 E_inf = -1.0 * cross(cfg_.u_inf, cfg_.B_imf);
    for (Id e = 0; e < cx.ent[1].n_owned; ++e)
    {
        if (bc_.edge_mode[e] == 1) r.E[e] = 0.0;                              // stationary perfectly conducting wall
        else if (bc_.edge_mode[e] == 2) r.E[e] = dot(E_inf, geo_.edge_vec[e]); // imposed free stream
    }
    emf_.finish(r.E.data());
}

void Solver::enforce_floors(State &s) const
{
    const CellComplex &cx = lm_.mesh.cx;
    for (int k = 0; k < 2; ++k)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            double *u = &s.U[k][5 * static_cast<std::size_t>(c)];
            const double ke = u[0] > 0 ? 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0] : 0.0;
            if (!(u[0] >= cfg_.rho_floor))
            {
                // density at/below the floor: restart that cell from a stagnant floor state (keeping the momentum would give u = m / floor)
                to_cons(Prim{cfg_.rho_floor, 0, 0, 0, cfg_.p_floor}, cfg_.sc.gamma, u);
            }
            else if ((cfg_.sc.gamma - 1.0) * (u[4] - ke) < cfg_.p_floor)
            {
                const Prim w = to_prim(u, cfg_.sc.gamma, cfg_.rho_floor, cfg_.p_floor);
                to_cons(w, cfg_.sc.gamma, u);
            }
        }
}

// out = base + dt * (rhs) on the owned fluid cells and the faces of owned cells; everything else is refreshed by the next halo exchange
void Solver::advance(const State &base, const Rhs &r, double dt, State &out)
{
    const CellComplex &cx = lm_.mesh.cx;
    out = base;
    for (int k = 0; k < 2; ++k)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            for (int m = 0; m < 5; ++m) out.U[k][5 * static_cast<std::size_t>(c) + m] += dt * r.dU[k][5 * static_cast<std::size_t>(c) + m];
        }
    faraday_update(cx, active_, r.E.data(), out.Phi.data(), dt);
    enforce_floors(out);
}

double Solver::step_first_order(double dt_max)
{
    const CellComplex &cx = lm_.mesh.cx;
    const double t_start = MPI_Wtime();
    prepare(s_);
    double gdt = local_dt();
    ex_.allreduce(&gdt, 1, MPI_MIN);
    gdt = std::min(gdt, dt_max);
    if (!std::isfinite(gdt) || gdt <= 0) throw std::runtime_error("invalid time step");
    const bool split = cfg_.sc.induce_coeff > 0 && cfg_.lorentz_force;
    if (split)
    {
        gyrate(lm_, cfg_, Btot_.data(), 0.5 * gdt, s_.U[0].data(), s_.U[1].data());
        prepare(s_);
    }
    dt_cur_ = gdt;
    rhs(s_, r0_); // fluid right-hand side and the explicit EMF (ideal + ambipolar, boundary values) at t^n
    // fluid: forward Euler
    for (int k = 0; k < 2; ++k)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            for (int m = 0; m < 5; ++m) s_.U[k][5 * static_cast<std::size_t>(c) + m] += gdt * r0_.dU[k][5 * static_cast<std::size_t>(c) + m];
        }
    enforce_floors(s_);
    // induction: implicit linear terms (resistivity, Hall) with the total EMF as the unknown, Phi only through d1(E)
    std::vector<double> E = r0_.E;
    if (ind_.active())
    {
        std::vector<double> ne(cx.n_cells(), 0.0);
        for (Id c = 0; c < cx.n_cells(); ++c)
            if (cx.cell_material[c] == kFluid) ne[c] = cfg_.sc.numinfo(W_[0][5 * static_cast<std::size_t>(c)], W_[1][5 * static_cast<std::size_t>(c)]).ne_eff;
        ind_.set_frozen_state(Btot_.data(), ne.data());
        std::vector<char> dir(cx.n_edges(), 0);
        for (Id e = 0; e < cx.n_edges(); ++e) dir[e] = bc_.edge_mode[e] == 2;
        double res = 0;
        cg_total_ += ind_.solve(s_.Phi.data(), r0_.E.data(), dir, gdt, E.data(), &res);
        ind_res_max_ = std::max(ind_res_max_, res);
    }
    faraday_update(cx, active_, E.data(), s_.Phi.data(), gdt);
    if (split) gyrate(lm_, cfg_, Btot_.data(), 0.5 * gdt, s_.U[0].data(), s_.U[1].data());
    time_ += gdt;
    ++nstep_;
    wall_s_ += MPI_Wtime() - t_start;
    return gdt;
}

double Solver::step(double dt_max)
{
    if (!initialized_) throw std::logic_error("Solver::step before initialize()");
    if (cfg_.time_order == 1) return step_first_order(dt_max);
    const double t_start = MPI_Wtime();
    prepare(s_);
    double gdt = local_dt();
    ex_.allreduce(&gdt, 1, MPI_MIN);
    gdt = std::min(gdt, dt_max);
    if (!std::isfinite(gdt) || gdt <= 0) throw std::runtime_error("invalid time step");
    // Strang splitting: exact half-step gyration, explicit RK2 for everything else, exact half-step gyration
    const bool split = cfg_.sc.induce_coeff > 0 && cfg_.lorentz_force;
    if (split)
    {
        gyrate(lm_, cfg_, Btot_.data(), 0.5 * gdt, s_.U[0].data(), s_.U[1].data());
        prepare(s_);
    }
    dt_cur_ = gdt;
    rhs(s_, r0_);
    advance(s_, r0_, gdt, s1_);
    prepare(s1_);
    rhs(s1_, r1_);
    advance(s1_, r1_, gdt, s2_);
    // SSP-RK2:  U^{n+1} = (U^n + U^{(2)}) / 2 ; the same combination of the CT updates keeps d2 B = 0
    for (int k = 0; k < 2; ++k)
        for (std::size_t i = 0; i < s_.U[k].size(); ++i) s_.U[k][i] = 0.5 * (s_.U[k][i] + s2_.U[k][i]);
    for (std::size_t i = 0; i < s_.Phi.size(); ++i) s_.Phi[i] = 0.5 * (s_.Phi[i] + s2_.Phi[i]);
    if (split) gyrate(lm_, cfg_, Btot_.data(), 0.5 * gdt, s_.U[0].data(), s_.U[1].data());
    // Mercury interior resistivity: implicit CT substep (operator splitting)
    if (res_.enabled()) cg_total_ += res_.advance(s_.Phi.data(), gdt, active_);
    time_ += gdt;
    ++nstep_;
    wall_s_ += MPI_Wtime() - t_start;
    return gdt;
}

Diagnostics Solver::collect(double dt)
{
    const CellComplex &cx = lm_.mesh.cx;
    prepare(s_);
    Diagnostics d;
    d.time = time_; d.step = nstep_; d.dt = dt;
    double loc[3] = {mercury::div_b(cx, s_.Phi.data(), PhiAdd_.data()).max_rel, 1e300, 1e300};
    double sum[2] = {0, 0};
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        loc[1] = std::min(loc[1], W_[0][5 * static_cast<std::size_t>(c)]);
        loc[2] = std::min(loc[2], W_[0][5 * static_cast<std::size_t>(c) + 4]);
        sum[0] += s_.U[0][5 * static_cast<std::size_t>(c)] * geo_.cell_volume[c];
        sum[1] += s_.U[1][5 * static_cast<std::size_t>(c)] * geo_.cell_volume[c];
    }
    double dth = 1e300;
    if (cfg_.sc.hall_coef > 0)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            const double r = norm(geo_.cell_center[c]);
            if (r < cfg_.hall_r_min) continue; // outside the radial taper
            double a = 0;
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k) a += geo_.face_area_mag[cx.c2f.idx[k]];
            const double h = 6.0 * geo_.cell_volume[c] / a;
            const double ne = cfg_.sc.numinfo(W_[0][5 * static_cast<std::size_t>(c)], W_[1][5 * static_cast<std::size_t>(c)]).ne_eff;
            const double B = norm(Vec3(Btot_[3 * c], Btot_[3 * c + 1], Btot_[3 * c + 2]));
            dth = std::min(dth, h * h * ne / (9.8696 * cfg_.sc.hall_coef * std::max(B, 1e-30)));
        }
    double g[3], gs[2];
    MPI_Allreduce(&loc[0], &g[0], 1, MPI_DOUBLE, MPI_MAX, comm_);
    MPI_Allreduce(&loc[1], &g[1], 2, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(sum, gs, 2, MPI_DOUBLE, MPI_SUM, comm_);
    MPI_Allreduce(&dth, &d.dt_hall, 1, MPI_DOUBLE, MPI_MIN, comm_);
    d.div_b = g[0]; d.min_rho = g[1]; d.min_p = g[2]; d.mass_H = gs[0]; d.mass_Na = gs[1];
    return d;
}

Diagnostics Solver::diagnostics() { return collect(0.0); }

Solver::Profile Solver::profile() const
{
    const CellComplex &cx = lm_.mesh.cx;
    Profile p;
    p.wall = wall_s_; p.comm = ex_.comm_seconds(); p.steps = nstep_; p.cg_iterations = cg_total_;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c) (cx.cell_material[c] == kSolid ? p.n_solid : p.n_fluid)++;
    p.n_ghost_cells = cx.n_cells() - cx.ent[3].n_owned;
    p.n_edges = cx.ent[1].n_owned; p.n_faces = cx.ent[2].n_owned; p.n_res_edges = res_.n_owned_edges();
    return p;
}
void Solver::reset_profile()
{
    wall_s_ = 0; cg_total_ = 0; ex_.reset_comm_clock();
}

void Solver::log(const Diagnostics &d) const
{
    if (rank_ == 0)
        std::printf("step %7ld  t=%.6e  dt=%.3e  divB=%.2e  min(rho_H)=%.3e  min(p_H)=%.3e  mass_H=%.8e  mass_Na=%.8e  dt_hall_explicit=%.2e\n", d.step, d.time, d.dt,
                    d.div_b, d.min_rho, d.min_p, d.mass_H, d.mass_Na, d.dt_hall);
}

void Solver::run()
{
    if (!initialized_) initialize();
    make_dir(cfg_.output_dir);
    log(collect(0.0));
    if (cfg_.write_vtk) write_vtk();
    double dt = 0;
    while (time_ < cfg_.t_end * (1.0 - 1e-12) && nstep_ < cfg_.max_steps)
    {
        dt = step(cfg_.t_end - time_);
        if (cfg_.log_every > 0 && nstep_ % cfg_.log_every == 0) log(collect(dt));
        if (cfg_.output_every > 0 && nstep_ % cfg_.output_every == 0 && nstep_ >= cfg_.output_start && cfg_.write_vtk) write_vtk();
        if (cfg_.restart_every > 0 && nstep_ % cfg_.restart_every == 0) write_global_fields(cfg_.output_dir + "/fields_final.unsb");
    }
    const Diagnostics d = collect(dt);
    log(d);
    if (cfg_.write_vtk) write_vtk();
    if (cfg_.write_fields) write_global_fields(cfg_.output_dir + "/fields_final.unsb");
    if (cfg_.write_profile) write_load_profile(cfg_.output_dir + "/load_profile");
}

// ------------------------------------------------------------------------------------------ output
void Solver::write_vtk()
{
    prepare(s_);
    const CellComplex &cx = lm_.mesh.cx;
    const Id n = cx.ent[3].n_owned;
    std::vector<io::VtkField> f = {{"rho_H", 1, {}}, {"rho_Na", 1, {}}, {"p_H", 1, {}}, {"p_Na", 1, {}}, {"u_H", 3, {}}, {"u_Na", 3, {}},
                                   {"B_induced", 3, {}}, {"B_total", 3, {}}, {"J", 3, {}}, {"u_plus", 3, {}}, {"material", 1, {}}};
    for (Id c = 0; c < n; ++c)
    {
        const std::size_t i = c;
        f[0].values.push_back(W_[0][5 * i]); f[1].values.push_back(W_[1][5 * i]);
        f[2].values.push_back(W_[0][5 * i + 4]); f[3].values.push_back(W_[1][5 * i + 4]);
        for (int a = 0; a < 3; ++a)
        {
            f[4].values.push_back(W_[0][5 * i + 1 + a]); f[5].values.push_back(W_[1][5 * i + 1 + a]);
            f[6].values.push_back(Bind_[3 * i + a]); f[7].values.push_back(Btot_[3 * i + a]);
            f[8].values.push_back(J_[3 * i + a]); f[9].values.push_back(uplus_[3 * i + a]);
        }
        f[10].values.push_back(cx.cell_material[c]);
    }
    char base[96];
    std::snprintf(base, sizeof base, "sol_%07ld", nstep_);
    char piece[160];
    std::snprintf(piece, sizeof piece, "%s/%s_r%05d.vtu", cfg_.output_dir.c_str(), base, rank_);
    make_dir(cfg_.output_dir);
    io::write_vtu(piece, lm_.mesh, n, f);
    MPI_Barrier(comm_);
    if (rank_ == 0)
    {
        std::vector<std::string> pieces;
        for (int r = 0; r < lm_.nranks; ++r)
        {
            char b[192];
            std::snprintf(b, sizeof b, "%s_r%05d.vtu", base, r);
            pieces.push_back(b);
        }
        std::vector<io::VtkField> proto;
        for (const auto &x : f) proto.push_back({x.name, x.ncomp, {}});
        io::write_pvtu(cfg_.output_dir + "/" + base + ".pvtu", pieces, proto);
        pvd_.emplace_back(time_, std::string(base) + ".pvtu");
        io::write_pvd(cfg_.output_dir + "/solution.pvd", pvd_);
    }
}

void Solver::write_global_fields(const std::string &path)
{
    prepare(s_);
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    std::vector<double> rho(nc), p(nc);
    for (std::size_t c = 0; c < nc; ++c) { rho[c] = W_[0][5 * c] + W_[1][5 * c]; p[c] = W_[0][5 * c + 4] + W_[1][5 * c + 4]; }
    struct Item { const char *name; Dim dim; int ncomp; const double *v; };
    const Item items[] = {{"U_H", Dim::Cell, 5, s_.U[0].data()}, {"U_Na", Dim::Cell, 5, s_.U[1].data()}, {"B", Dim::Face, 1, s_.Phi.data()},
                          {"rho", Dim::Cell, 1, rho.data()},      {"p", Dim::Cell, 1, p.data()},        {"B_cell", Dim::Cell, 3, Btot_.data()},
                          {"J_cell", Dim::Cell, 3, J_.data()}};
    std::vector<io::FieldRecord> out;
    for (const Item &it : items)
    {
        auto g = gather_entities(comm_, cx.ent[di(it.dim)], it.v, it.ncomp);
        if (rank_ == 0) out.push_back({it.name, it.dim, it.ncomp, std::move(g)});
    }
    if (rank_ == 0)
    {
        out.push_back({"_time", Dim::Node, 1, {time_}});
        // hash of the global complex is not available on a sub-complex run; 0 marks "not checked"
        io::write_fields(path, 0, out);
        std::printf("wrote %s (t = %.6e)\n", path.c_str(), time_);
    }
}
} // namespace mercury

// ------------------------------------------------------------------------------------------ load-balance profile
namespace mercury
{
void Solver::write_load_profile(const std::string &base)
{
    const Profile p = profile();
    const int np = lm_.nranks;
    const double v[9] = {p.wall, p.comm, double(p.n_fluid), double(p.n_solid), double(p.n_ghost_cells), double(p.n_res_edges), double(p.cg_iterations),
                         double(p.n_edges), double(p.steps)};
    std::vector<double> all(9 * static_cast<std::size_t>(np));
    MPI_Gather(v, 9, MPI_DOUBLE, all.data(), 9, MPI_DOUBLE, 0, comm_);
    if (rank_ != 0) return;
    auto at = [&](int r, int k) { return all[9 * static_cast<std::size_t>(r) + k]; };
    const double steps = std::max(at(0, 8), 1.0);
    // per-step compute time of each rank = wall - communication/wait
    std::vector<double> tc(np), tw(np);
    double sf = 0, ss = 0, sff = 0, sfs = 0, sss = 0, stf = 0, sts = 0, tmax = 0, tsum = 0;
    for (int r = 0; r < np; ++r)
    {
        tc[r] = (at(r, 0) - at(r, 1)) / steps;
        tw[r] = at(r, 1) / steps;
        const double nf = at(r, 2), ns = at(r, 3);
        sff += nf * nf; sfs += nf * ns; sss += ns * ns; stf += nf * tc[r]; sts += ns * tc[r];
        sf += nf; ss += ns; tmax = std::max(tmax, tc[r]); tsum += tc[r];
    }
    // least squares  tc ~ a nf + b ns   (non-negative)
    double a = 0, b = 0;
    const double det = sff * sss - sfs * sfs;
    if (det > 1e-9 * sff * sss) { a = (stf * sss - sts * sfs) / det; b = (sts * sff - stf * sfs) / det; }
    bool fitted = det > 1e-9 * sff * sss && a > 0 && b >= 0;
    if (!fitted) { a = sff > 0 ? stf / sff : 0; b = 0; fitted = false; } // fluid-only fit; solid ratio kept at the default
    double res = 0, tot = 0;
    for (int r = 0; r < np; ++r)
    {
        const double pred = a * at(r, 2) + b * at(r, 3);
        res += (tc[r] - pred) * (tc[r] - pred);
        tot += (tc[r] - tsum / np) * (tc[r] - tsum / np);
    }
    const double r2 = tot > 0 ? 1.0 - res / tot : 0.0;
    // communication/wait: seconds per step per ghost cell (median over ranks), in units of one fluid-cell update
    std::vector<double> cg;
    for (int r = 0; r < np; ++r) if (at(r, 4) > 0) cg.push_back(tw[r] / at(r, 4));
    std::sort(cg.begin(), cg.end());
    const double wait_per_ghost = cg.empty() ? 0.0 : cg[cg.size() / 2];
    CostModelRecommendation rec;
    rec.cost_fluid = 1.0;
    rec.cost_solid = (fitted && a > 0) ? std::min(5.0, std::max(0.05, b / a)) : 0.35;
    rec.comm_fluid = a > 0 ? std::min(5.0, wait_per_ghost / a) : 1.0;
    rec.comm_solid = 0.3 * rec.comm_fluid;
    const double mean = tsum / np;
    double wsum = 0, tw_all = 0;
    for (int r = 0; r < np; ++r) { wsum += at(r, 1); tw_all += at(r, 0); }
    const double wait_frac = tw_all > 0 ? wsum / tw_all : 0.0;
    FILE *fj = std::fopen((base + ".json").c_str(), "w");
    if (!fj) return;
    std::fprintf(fj, "{\n  \"description\": \"MercuryUns load profile of a coarse run; cost_* feed mpcns_pre (cost_profile=)\",\n");
    std::fprintf(fj, "  \"nranks\": %d, \"steps\": %.0f,\n  \"fit\": {\"seconds_per_fluid_cell_per_step\": %.4e, \"seconds_per_solid_cell_per_step\": %.4e, "
                     "\"solid_to_fluid_ratio\": %.4f, \"two_parameter_fit_ok\": %s, \"r2\": %.4f, \"wait_seconds_per_ghost_cell_per_step\": %.4e},\n",
                 np, steps, a, b, a > 0 ? b / a : 0.0, fitted ? "true" : "false", r2, wait_per_ghost);
    std::fprintf(fj, "  \"resistive_cg_iterations_per_step\": %.2f,\n", at(0, 6) / steps);
    std::fprintf(fj, "  \"measured_imbalance\": {\"compute_max_over_mean\": %.4f, \"wait_fraction_of_wall\": %.4f},\n", mean > 0 ? tmax / mean : 1.0,
                 [&] { double w = 0, t = 0; for (int r = 0; r < np; ++r) { w += at(r, 1); t += at(r, 0); } return t > 0 ? w / t : 0.0; }());
    std::fprintf(fj, "  \"recommended_cost_model\": {\"cost_fluid\": %.4f, \"cost_solid\": %.4f, \"comm_fluid\": %.4f, \"comm_solid\": %.4f},\n", rec.cost_fluid,
                 rec.cost_solid, rec.comm_fluid, rec.comm_solid);
    std::fprintf(fj, "  \"ranks\": [\n");
    for (int r = 0; r < np; ++r)
        std::fprintf(fj, "    {\"rank\": %d, \"compute_s_per_step\": %.4e, \"wait_s_per_step\": %.4e, \"wait_fraction\": %.3f, \"fluid_cells\": %.0f, \"solid_cells\": %.0f, "
                         "\"ghost_cells\": %.0f, \"ghost_to_owned\": %.3f, \"owned_edges\": %.0f, \"resistive_edges\": %.0f, \"cg_iterations\": %.0f}%s\n",
                     r, tc[r], tw[r], at(r, 0) > 0 ? at(r, 1) / at(r, 0) : 0.0, at(r, 2), at(r, 3), at(r, 4), (at(r, 2) + at(r, 3)) > 0 ? at(r, 4) / (at(r, 2) + at(r, 3)) : 0.0,
                     at(r, 7), at(r, 5), at(r, 6), r + 1 < np ? "," : "");
    std::fprintf(fj, "  ]\n}\n");
    std::fclose(fj);
    FILE *ft = std::fopen((base + ".txt").c_str(), "w");
    if (ft)
    {
        std::fprintf(ft, "// written by MercuryUns (%d ranks, %.0f steps): cost model for mpcns_pre, use  cost_profile=%s.txt\n", np, steps, base.c_str());
        std::fprintf(ft, "double cost_fluid %.4f\ndouble cost_solid %.4f\ndouble comm_fluid %.4f\ndouble comm_solid %.4f\n", rec.cost_fluid, rec.cost_solid, rec.comm_fluid, rec.comm_solid);
        std::fclose(ft);
    }
    std::printf("load profile: compute max/mean = %.3f, wait fraction = %.3f, solid/fluid cost = %.3f (%s fit, r2 = %.2f), comm per ghost cell = %.3f cells -> %s.json/.txt\n",
                mean > 0 ? tmax / mean : 1.0, wait_frac, rec.cost_solid, fitted ? "2-parameter" : "default-ratio", r2, rec.comm_fluid, base.c_str());
}
} // namespace mercury
