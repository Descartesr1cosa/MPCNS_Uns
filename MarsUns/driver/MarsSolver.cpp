#include "driver/MarsSolver.h"

#include <filesystem>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "diagnostics/DivB.h"
#include "induction/Current.h"
#include "induction/Faraday.h"
#include "io/GlobalFields.h"
#include "physics/Euler.h"
#include "uns/geometry/Reconstruct.h"
#include "uns/io/MeshIO.h"
#include "uns/io/VtkWriter.h"

namespace mars
{
using namespace mhd;
using namespace uns;

namespace
{
void make_dir(const std::string &d)
{
    std::error_code ec; // mkdir -p, portable
    std::filesystem::create_directories(d, ec);
}
inline std::size_t S(Id c) { return static_cast<std::size_t>(c); }
const char *name_of(int s) { return MARS::SpeciesNames[s].data(); }
} // namespace

// ------------------------------------------------------------------------------------------------ parameters
MarsParams MarsParams::from_params(const Params &p)
{
    MarsParams m;
    m.cfg = Config::from_params(p);
    const double gamma = p.get_double("gamma", 1.4);
    m.NA = p.get_double("NA", 6.02214076e23);
    const double R = p.get_double("R_uni", 8.314462618);
    m.mu0 = p.get_double("mu_mag", 1.256637061e-6);
    m.q = p.get_double("q_e", 1.602176634e-19);
    m.kB = R / m.NA;
    m.U = p.get_double("U_ref", 3.27e4);
    m.T = p.get_double("T_ref", 1.286e5);
    m.n = p.get_double("n_ref", 4.0e6);
    m.B = p.get_double("B_ref", 3.0e-9);
    m.L = p.get_double("L_ref", 3.396e6);
    const double Mref = p.get_double("Molecular_mass", 0.001);
    m.rho_ref = Mref / m.NA * m.n;
    m.p_ref = m.rho_ref * m.U * m.U;
    m.time_ref = m.L / m.U;
    m.hall_scale = m.B / (m.mu0 * m.L * m.q * m.n * m.U);
    m.pressure_scale = m.p_ref / (m.q * m.n * m.L * m.U * m.B);
    m.mag_scale = m.B * m.B / (m.mu0 * m.rho_ref * m.U * m.U);
    m.ne_floor = p.get_double("ne_hall_floor", 1e-2);
    m.rho_floor = p.get_double("rho_floor", 1e-10);
    m.p_floor = p.get_double("p_floor", 1e-10);

    // solar wind (the INITIAL / Boundary_Farfield lists of the original setup): physical values
    const double sw_n = p.get_double("sw_n", m.n), sw_T = p.get_double("sw_T", 1.75e5), sw_U = p.get_double("sw_U", 5.0e5);
    const double sw_M = p.get_double("sw_Molecular_mass", Mref);
    const double cy = p.get_double("sw_cy", 0.0), cz = p.get_double("sw_cz", 0.0);
    const double cx = -std::sqrt(std::max(0.0, 1.0 - cy * cy - cz * cz));
    m.sw_rho = (sw_n / m.n) * (sw_M / Mref);
    m.sw_u = {cx * sw_U / m.U, cy * sw_U / m.U, cz * sw_U / m.U};
    m.sw_p = sw_n * m.kB * sw_T / m.p_ref;
    m.sw_B = {p.get_double("sw_Bx", -1.68e-9) / m.B, p.get_double("sw_By", 2.49e-9) / m.B, p.get_double("sw_Bz", 0.0) / m.B};
    m.wall_T = p.get_double("wall_T", 1500.0);
    m.hall = p.get_bool("is_Hall_Efield", true);
    m.ambipolar = p.get_bool("is_Ambipolar_Efield", true);
    m.crustal = p.get_bool("is_Crustal_Bfield", false);
    m.dipole_file = p.get_string("Crustal_Dipoles_File", "Dipoles");
    m.cutoff_km = p.get_double("Crustal_Dipoles_Cutoff_km", 0.0);
    m.atm.mars_radius_m = p.get_double("mars_radius_m", m.L);
    m.atm.lower_boundary_altitude_km = p.get_double("mars_lower_boundary_km", 100.0);
    m.atm.heliocentric_distance_au = p.get_double("mars_orbit_au", 1.524);
    m.atm.sun_direction = {p.get_double("mars_sun_x", 1.0), p.get_double("mars_sun_y", 0.0), p.get_double("mars_sun_z", 0.0)};
    m.atm.ion_temperature_K = m.atm.electron_temperature_K = m.wall_T;

    // the generic configuration used by the shared Mercury modules
    Config &c = m.cfg;
    c.sc.gamma = gamma;
    c.sc.hall_coef = m.hall_scale;
    c.sc.pressure_coef = m.pressure_scale;
    c.sc.inver_MA2 = m.mag_scale;
    c.sc.ne_floor = m.ne_floor;
    c.cfl = p.get_double("CFL", 0.5);
    c.second_order = p.get_bool("second_order", false); // MarsZ5 is first order
    c.rho_floor = m.rho_floor;
    c.p_floor = m.p_floor;
    c.u_inf = m.sw_u;
    c.B_imf = m.sw_B;
    c.hall = m.hall;
    c.ambipolar = m.ambipolar;
    c.hall_r_min = c.hall_r_max = 0.0; // no planetary taper: the Hall term acts everywhere
    c.resistive = false;               // no Mercury interior resistivity
    c.chemistry = false;               // the Mercury Na+ source terms are not part of this branch
    c.time_order = 1;
    c.steady = false;
    c.implicit_ideal = false;
    c.mg_levels = 0;
    c.dipole = false;
    return m;
}

// ------------------------------------------------------------------------------------------------ construction
MarsSolver::MarsSolver(MPI_Comm comm, const LocalMesh &lm, const MarsParams &prm)
    : comm_(comm), lm_(lm), prm_(prm), cfg_(prm.cfg), geo_(compute_geometry(lm.mesh)), ex_(comm, lm), emf_asm_(lm, geo_, ex_),
      bc_(BoundaryInfo::build(lm, geo_, cfg_)), stencil_(lm, geo_, bc_), ind_(lm, geo_, ex_, cfg_), active_(active_faces(lm.mesh.cx)),
      atmosphere_(prm.atm)
{
    MPI_Comm_rank(comm, &rank_);
    stencil_.set_rel_limit(cfg_.recon_rel_limit);
    for (int s = 0; s < NS; ++s)
        fs_[s] = (s == MARS::H) ? Prim{prm_.sw_rho, prm_.sw_u.x, prm_.sw_u.y, prm_.sw_u.z, prm_.sw_p}
                                : Prim{prm_.rho_floor, 0, 0, 0, std::max(prm_.p_floor, prm_.sw_p * 1e-6)};
    ghost_ = [this](Id f, int s, const Prim &in, const Vec3 &n) -> Prim {
        const FaceKind k = bc_.face_kind[f];
        if (k == FaceKind::Wall) return wall_state_(f, s);
        return ghost_state(k, s, in, n, fs_[s]);
    };
    // node -> fluid cells (nodal electron pressure of the ambipolar EMF)
    const CellComplex &cx = lm.mesh.cx;
    nptr_.assign(cx.n_nodes() + 1, 0);
    std::vector<std::pair<Id, Id>> pr;
    std::vector<Id> cn;
    for (Id c = 0; c < cx.n_cells(); ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        cx.cell_nodes(c, cn);
        for (Id v : cn) pr.emplace_back(v, c);
    }
    for (auto &q : pr) ++nptr_[q.first + 1];
    for (Id v = 0; v < cx.n_nodes(); ++v) nptr_[v + 1] += nptr_[v];
    ncell_.resize(pr.size());
    std::vector<std::int64_t> fill(nptr_.begin(), nptr_.end() - 1);
    for (auto &q : pr) ncell_[fill[q.first]++] = q.second;
}

// Dirichlet ghost state of the inner (photochemical-equilibrium) boundary, evaluated at the boundary face centre.
Prim MarsSolver::wall_state_(Id f, int s)
{
    auto it = wall_cache_.find(f);
    if (it == wall_cache_.end())
    {
        const Vec3 x = prm_.L * geo_.face_center[f];
        const MARS::PceSolution sol = atmosphere_.WallIonDensity({x.x, x.y, x.z});
        std::array<Prim, NS> st;
        for (int k = 0; k < NS; ++k)
        {
            double rho = sol.ion_cm3[k] * 1.0e6 * MARS::ParticleMass[k] / prm_.rho_ref;
            if (k == MARS::H) rho = 0.3 * prm_.sw_rho; // the original H+ condition: 0.3 x upstream density
            rho = std::max(rho, prm_.rho_floor);
            const double number = rho * prm_.rho_ref / MARS::ParticleMass[k];
            st[k] = Prim{rho, 0, 0, 0, std::max(number * prm_.kB * prm_.wall_T / prm_.p_ref, prm_.p_floor)};
        }
        it = wall_cache_.emplace(f, st).first;
    }
    return it->second[s];
}

void MarsSolver::build_background_()
{
    const CellComplex &cx = lm_.mesh.cx;
    PhiAdd_.assign(cx.n_faces(), 0.0);
    BaddCell_.assign(3 * S(cx.n_cells()), 0.0);
    for (Id f = 0; f < cx.n_faces(); ++f) PhiAdd_[f] = dot(prm_.sw_B, geo_.face_area[f]); // uniform IMF
    if (prm_.crustal)
    {
        std::ifstream in(prm_.dipole_file);
        if (!in) throw std::runtime_error("cannot open crustal dipole file " + prm_.dipole_file);
        MARS::DipoleModel model;
        const std::string src((std::istreambuf_iterator<char>(in)), {});
        std::istringstream rec(src);
        model.Read(rec, prm_.L, prm_.B, prm_.mu0, prm_.dipole_file);
        model.SetCutoff(prm_.cutoff_km * 1000.0 / prm_.L);
        // A.dl along every edge (identical on all copies of an edge: it depends only on its two node coordinates), then B = d1 A
        std::vector<double> A(cx.n_edges(), 0.0);
        for (Id e = 0; e < cx.n_edges(); ++e)
        {
            Id a, b;
            cx.edge_nodes(e, a, b);
            const Vec3 &xa = lm_.mesh.x[a], &xb = lm_.mesh.x[b];
            A[e] = model.EdgeIntegral({xa.x, xa.y, xa.z}, {xb.x, xb.y, xb.z});
        }
        for (Id f = 0; f < cx.n_faces(); ++f)
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) PhiAdd_[f] += cx.f2e.sgn[k] * A[cx.f2e.idx[k]];
        if (rank_ == 0) std::printf("[Mars Dipoles] %zu sources from %s, cutoff %.0f km\n", model.size(), prm_.dipole_file.c_str(), prm_.cutoff_km);
    }
    else if (rank_ == 0)
        std::printf("[Mars Dipoles] disabled; the IMF is still part of the background field\n");
    const auto b = reconstruct_cell_vector(cx, geo_, PhiAdd_.data(), cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c) { BaddCell_[3 * S(c)] = b[c].x; BaddCell_[3 * S(c) + 1] = b[c].y; BaddCell_[3 * S(c) + 2] = b[c].z; }
}

void MarsSolver::build_atmosphere_()
{
    const CellComplex &cx = lm_.mesh.cx;
    neutral_.assign(S(cx.n_cells()) * MARS::NumNeutral, 0.0);
    photo_.assign(S(cx.n_cells()) * MARS::NumPhoto, 0.0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const Vec3 x = prm_.L * geo_.cell_center[c];
        const MARS::AtmosphereState st = atmosphere_.Evaluate({x.x, x.y, x.z});
        for (std::size_t n = 0; n < MARS::NumNeutral; ++n) neutral_[S(c) * MARS::NumNeutral + n] = st.neutral_cm3[n];
        for (std::size_t n = 0; n < MARS::NumPhoto; ++n) photo_[S(c) * MARS::NumPhoto + n] = st.photo_cm3_s[n];
    }
}

void MarsSolver::initialize()
{
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    for (int s = 0; s < NS; ++s)
    {
        U_[s].assign(5 * nc, 0.0);
        W_[s].assign(5 * nc, 0.0);
        grad_[s].assign(15 * nc, 0.0);
        rhs_[s].assign(5 * nc, 0.0);
        for (std::size_t c = 0; c < nc; ++c) to_cons(fs_[s], cfg_.sc.gamma, &U_[s][5 * c]);
    }
    Phi_.assign(cx.n_faces(), 0.0);
    if (!cfg_.init_fields.empty())
    {
        const auto fields = io::read_fields(cfg_.init_fields);
        int got = 0;
        bool gotB = false;
        for (const auto &fr : fields)
        {
            if (fr.name == "_time" && !fr.values.empty()) { time_ = fr.values[0]; continue; }
            if (fr.name == "B" && fr.dim == Dim::Face)
            {
                for (Id f = 0; f < cx.n_faces(); ++f) Phi_[f] = fr.values.at(S(cx.ent[2].gid[f]));
                gotB = true;
                continue;
            }
            for (int s = 0; s < NS; ++s)
                if (fr.name == std::string("U_") + name_of(s) && fr.dim == Dim::Cell && fr.ncomp == 5)
                {
                    for (std::size_t c = 0; c < nc; ++c)
                        for (int m = 0; m < 5; ++m) U_[s][5 * c + m] = fr.values.at(S(cx.ent[3].gid[c]) * 5 + m);
                    ++got;
                }
        }
        if (got != NS || !gotB) throw std::runtime_error("init_fields '" + cfg_.init_fields + "' lacks one of the ten U_<ion> fields or B");
    }
    PhiTot_.assign(cx.n_faces(), 0.0);
    Bind_.assign(3 * nc, 0.0);
    Btot_.assign(3 * nc, 0.0);
    J_.assign(3 * nc, 0.0);
    uplus_.assign(3 * nc, 0.0);
    dpe_.assign(3 * nc, 0.0);
    ne_.assign(nc, 0.0);
    ne_eff_.assign(nc, 0.0);
    pe_.assign(nc, 0.0);
    build_background_();
    build_atmosphere_();
    initialized_ = true;
    prepare();
}

// ------------------------------------------------------------------------------------------------ derived fields
void MarsSolver::prepare()
{
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    for (int s = 0; s < NS; ++s) ex_.forward(Dim::Cell, U_[s].data(), 5);
    ex_.forward(Dim::Face, Phi_.data(), 1);
    const double g = cfg_.sc.gamma;
    for (int s = 0; s < NS; ++s)
        for (std::size_t c = 0; c < nc; ++c)
        {
            const Prim w = to_prim(&U_[s][5 * c], g, cfg_.rho_floor, cfg_.p_floor);
            double *o = &W_[s][5 * c];
            o[0] = w.rho; o[1] = w.u; o[2] = w.v; o[3] = w.w; o[4] = w.p;
        }
    const auto bi = reconstruct_cell_vector(cx, geo_, Phi_.data(), cx.n_cells());
    for (std::size_t c = 0; c < nc; ++c)
    {
        const double b[3] = {bi[c].x, bi[c].y, bi[c].z};
        for (int a = 0; a < 3; ++a) { Bind_[3 * c + a] = b[a]; Btot_[3 * c + a] = b[a] + BaddCell_[3 * c + a]; }
        double ne = 0, pe = 0, v[3] = {0, 0, 0};
        for (int s = 0; s < NS; ++s)
        {
            const double n = std::max(U_[s][5 * c], 0.0) * prm_.rho_ref / (MARS::ParticleMass[s] * prm_.n);
            ne += n;
            pe += W_[s][5 * c + 4];
            for (int a = 0; a < 3; ++a) v[a] += n * W_[s][5 * c + 1 + a];
        }
        ne_[c] = ne;
        ne_eff_[c] = std::sqrt(ne * ne + prm_.ne_floor * prm_.ne_floor);
        pe_[c] = pe;
        const double inv = ne > 1e-300 ? 1.0 / ne : 0.0;
        for (int a = 0; a < 3; ++a) uplus_[3 * c + a] = (cx.cell_material[c] == kFluid) ? v[a] * inv : 0.0;
    }
    for (Id f = 0; f < cx.n_faces(); ++f) PhiTot_[f] = Phi_[f] + PhiAdd_[f];
    cell_current(lm_, geo_, Bind_.data(), J_.data());
    // grad p_e by Gauss' theorem with face-averaged pressure
    std::fill(dpe_.begin(), dpe_.end(), 0.0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        Vec3 sum;
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            double pf = pe_[c];
            if (cx.f2c.count(f) == 2)
            {
                const Id o = cx.f2c.idx[cx.f2c.begin(f)] == c ? cx.f2c.idx[cx.f2c.begin(f) + 1] : cx.f2c.idx[cx.f2c.begin(f)];
                if (cx.cell_material[o] == kFluid) pf = 0.5 * (pf + pe_[o]);
            }
            sum += (pf * static_cast<double>(cx.c2f.sgn[k])) * geo_.face_area[f];
        }
        const double iv = 1.0 / geo_.cell_volume[c];
        dpe_[3 * S(c)] = sum.x * iv; dpe_[3 * S(c) + 1] = sum.y * iv; dpe_[3 * S(c) + 2] = sum.z * iv;
    }
}

// ------------------------------------------------------------------------------------------------ time step
double MarsSolver::local_dt_() const
{
    const CellComplex &cx = lm_.mesh.cx;
    double dt = 1e300;
    lim_ = DtLimit();
    auto note = [&](double d, Id c, int kind) { if (d < lim_.dt) { lim_.dt = d; lim_.x = geo_.cell_center[c]; lim_.kind = kind; } };
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const double V = geo_.cell_volume[c];
        // (1) each ion on its own: cfl V / sum_f (|u.A| + c_s |A|)
        for (int s = 0; s < NS; ++s)
        {
            const double *w = &W_[s][5 * S(c)];
            const double cs = std::sqrt(cfg_.sc.gamma * w[4] / w[0]);
            double den = 0;
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
            {
                const Id f = cx.c2f.idx[k];
                const Vec3 &A = geo_.face_area[f];
                den += std::abs(w[1] * A.x + w[2] * A.y + w[3] * A.z) + cs * geo_.face_area_mag[f];
            }
            if (den > 0) { dt = std::min(dt, cfg_.cfl * V / den); note(cfg_.cfl * V / den, c, s); }
        }
        // (2) the charged plasma as one fluid: fast magnetosonic waves of the total field, u+ as the advection velocity
        double rho = 0, p = 0;
        for (int s = 0; s < NS; ++s) { rho += std::max(U_[s][5 * S(c)], 0.0); p += std::max(W_[s][5 * S(c) + 4], 0.0); }
        rho = std::max(rho, prm_.rho_floor);
        p = std::max(p, prm_.p_floor);
        const Vec3 b(Btot_[3 * S(c)], Btot_[3 * S(c) + 1], Btot_[3 * S(c) + 2]);
        const Vec3 up(uplus_[3 * S(c)], uplus_[3 * S(c) + 1], uplus_[3 * S(c) + 2]);
        const double cs2 = cfg_.sc.gamma * p / rho, va2 = prm_.mag_scale * dot(b, b) / rho;
        double den = 0;
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            const double A = geo_.face_area_mag[f];
            const Vec3 &Av = geo_.face_area[f];
            const double bn = dot(b, Av) / A, van2 = prm_.mag_scale * bn * bn / rho;
            const double sum = cs2 + va2, disc = std::max(0.0, sum * sum - 4.0 * cs2 * van2);
            const double cf = std::sqrt(0.5 * (sum + std::sqrt(disc)));
            den += std::abs(dot(up, Av)) + cf * A;
        }
        if (den > 0) { dt = std::min(dt, cfg_.cfl * V / den); note(cfg_.cfl * V / den, c, NS); }
    }
    return dt;
}

double MarsSolver::time_step_limit()
{
    prepare();
    double d = local_dt_();
    ex_.allreduce(&d, 1, MPI_MIN);
    return d;
}

// ------------------------------------------------------------------------------------------------ fluxes and sources
void MarsSolver::fluid_and_sources_(double)
{
    const CellComplex &cx = lm_.mesh.cx;
    for (int s = 0; s < NS; ++s)
    {
        if (cfg_.second_order)
        {
            stencil_.limited_gradients(lm_, geo_, bc_, W_[s].data(), s, fs_[s], true, grad_[s].data(), &ghost_);
            ex_.forward(Dim::Cell, grad_[s].data(), 15);
        }
        else
            std::fill(grad_[s].begin(), grad_[s].end(), 0.0);
        fluid_flux_rhs(lm_, geo_, bc_, cfg_, s, fs_[s], W_[s].data(), grad_[s].data(), rhs_[s].data(), &ghost_);
    }
    // Explicit part of the MarsZ5 sources: J x B and the ambipolar force (weights w_s = n_s / n_e,eff).  The stiff parts -- ion gyration
    // about B, ion-neutral drag with its frictional heating, thermal relaxation to the neutrals -- are point-implicit, see implicit_sources_().
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const std::size_t i = c;
        const Vec3 b(Btot_[3 * i], Btot_[3 * i + 1], Btot_[3 * i + 2]), j(J_[3 * i], J_[3 * i + 1], J_[3 * i + 2]);
        const Vec3 dpe(dpe_[3 * i], dpe_[3 * i + 1], dpe_[3 * i + 2]);
        const Vec3 jxb = cross(j, b);
        for (int s = 0; s < NS; ++s)
        {
            const double rho = std::max(U_[s][5 * i], 0.0);
            const double w = rho * prm_.rho_ref / (MARS::ParticleMass[s] * prm_.n) / ne_eff_[i]; // charge share  n_s / n_e,eff
            const double *wp = &W_[s][5 * i];
            Vec3 f = (w * prm_.mag_scale) * jxb;
            if (cfg_.ambipolar) f -= w * dpe;
            double *r = &rhs_[s][5 * i];
            r[1] += f.x; r[2] += f.y; r[3] += f.z;
            r[4] += dot(f, Vec3(wp[1], wp[2], wp[3]));
        }
    }
}

// Point-implicit (cell-local) treatment of the stiff sources, applied to U after the explicit update:
//  1. gyration  dz_s/dt = -i Omega_s (z_s - Z),  z = u_x + i u_y in the plane normal to B, Z = z+ = sum chi_s z_s (u_par is unchanged).
//     Crank-Nicolson, solved in closed form for the coupled ten-ion system:  Z' = P / (1 - Q), P and Q below.  A-stable, second order, conserves
//     the mass momentum sum rho_s u_s exactly (Z is the number-weighted mean, so sum n_s (z_s - Z) = 0 at both time levels); the kinetic-energy change is added to E (the thermal energy is untouched).
//  2. drag  dm/dt = -nu m, nu = t_ref sum_n n_n 2e-10  (backward Euler);  the lost kinetic energy stays in E, i.e. becomes heat (as in MarsZ5).
//  3. thermal relaxation + frictional heating  dp/dt = (gamma-1) [k (185 - T) + S_f],  k = a5 n_s sum_n m_s nu_n/(m_s+m_n)  (backward Euler in p).
// None of these limits the time step any more.
void MarsSolver::implicit_sources_(double dt)
{
    const CellComplex &cx = lm_.mesh.cx;
    const double g = cfg_.sc.gamma;
    const double a5 = 3.0e6 * prm_.L * prm_.kB / (prm_.rho_ref * prm_.U * prm_.U * prm_.U);
    using cd = std::complex<double>;
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const std::size_t i = c;
        double nrel[NS], ntot = 0;
        for (int s = 0; s < NS; ++s)
        {
            nrel[s] = std::max(U_[s][5 * i], 0.0) * prm_.rho_ref / (MARS::ParticleMass[s] * prm_.n);
            ntot += nrel[s];
        }
        // ---- 1. gyration
        const Vec3 B(Btot_[3 * i], Btot_[3 * i + 1], Btot_[3 * i + 2]);
        const double Bm = norm(B);
        if (ntot > 1e-300 && Bm > 1e-30)
        {
            const Vec3 k = (1.0 / Bm) * B;
            const Vec3 t = std::abs(k.x) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
            Vec3 e1 = cross(k, t);
            e1 = (1.0 / norm(e1)) * e1;
            const Vec3 e2 = cross(k, e1);
            cd z[NS], al[NS], Z0 = 0;
            for (int s = 0; s < NS; ++s)
            {
                const double *u = &U_[s][5 * i];
                const double rho = std::max(u[0], prm_.rho_floor);
                const Vec3 v(u[1] / rho, u[2] / rho, u[3] / rho);
                z[s] = cd(dot(v, e1), dot(v, e2));
                al[s] = cd(0.0, prm_.q * prm_.B * prm_.L / (MARS::ParticleMass[s] * prm_.U) * Bm * dt * 0.5); // i alpha_s / 2
                Z0 += (nrel[s] / ntot) * z[s];
            }
            // dz/dt = -i Omega (z - Z): the CN coefficient is a = -i alpha / 2 (al holds +i alpha / 2)
            cd Pn = 0, Qn = 0;
            for (int s = 0; s < NS; ++s)
            {
                const double chi = nrel[s] / ntot;
                const cd a = -al[s];
                Pn += chi * (z[s] * (1.0 - a) + a * Z0) / (1.0 + a);
                Qn += chi * a / (1.0 + a);
            }
            const cd Zn = Pn / (1.0 - Qn);
            for (int s = 0; s < NS; ++s)
            {
                const cd a = -al[s];
                const cd zn = (z[s] * (1.0 - a) + a * (Z0 + Zn)) / (1.0 + a);
                double *u = &U_[s][5 * i];
                const double rho = std::max(u[0], prm_.rho_floor);
                const Vec3 v(u[1] / rho, u[2] / rho, u[3] / rho);
                const double vpar = dot(v, k);
                const Vec3 vn = vpar * k + zn.real() * e1 + zn.imag() * e2;
                u[4] += 0.5 * rho * (dot(vn, vn) - dot(v, v));
                u[1] = rho * vn.x; u[2] = rho * vn.y; u[3] = rho * vn.z;
            }
        }
        // ---- 2./3. drag and thermal relaxation
        double coll_n[MARS::NumReactiveNeutral], collision = 0;
        for (std::size_t n = 0; n < MARS::NumReactiveNeutral; ++n)
        {
            coll_n[n] = std::max(neutral_[i * MARS::NumNeutral + n], 0.0) * 2.0e-10;
            collision += coll_n[n];
        }
        const double nu = prm_.time_ref * collision;
        for (int s = 0; s < NS; ++s)
        {
            double *u = &U_[s][5 * i];
            const double m = MARS::ParticleMass[s];
            const double rho = std::max(u[0], prm_.rho_floor);
            const Vec3 v(u[1] / rho, u[2] / rho, u[3] / rho);
            const double ke_old = 0.5 * rho * dot(v, v);
            const double den = 1.0 + nu * dt;
            u[1] /= den; u[2] /= den; u[3] /= den;
            const double ke = ke_old / (den * den);
            double p = std::max((g - 1.0) * (u[4] - ke), prm_.p_floor); // thermal pressure incl. the heat of the drag
            double b1c = 0, b2 = 0;
            for (std::size_t n = 0; n < MARS::NumReactiveNeutral; ++n)
            {
                const double mm = m + MARS::ParticleMass[n];
                b1c += m * coll_n[n] / mm;
                b2 += MARS::ParticleMass[n] * coll_n[n] / mm;
            }
            const double n_m3 = rho * prm_.rho_ref / m;
            const double ns_cm3 = n_m3 * 1e-6;
            const double kk = a5 * ns_cm3 * b1c;
            const double cT = prm_.p_ref / (std::max(n_m3, 1.0) * prm_.kB);
            const double Sf = prm_.time_ref * rho * dot(v, v) * b2;
            double pn;
            if (p * cT >= 50.0) pn = (p + dt * (g - 1.0) * (kk * 185.0 + Sf)) / (1.0 + dt * (g - 1.0) * kk * cT);
            else pn = p + dt * (g - 1.0) * (kk * (185.0 - 50.0) + Sf);
            u[4] = pn / (g - 1.0) + ke;
        }
    }
}

// ------------------------------------------------------------------------------------------------ electromotive force
void MarsSolver::ambipolar_emf_(double *E) const
{
    const CellComplex &cx = lm_.mesh.cx;
    std::vector<double> pn(cx.n_nodes(), 0.0);
    std::vector<char> have(cx.n_nodes(), 0);
    auto node_p = [&](Id v) {
        if (have[v]) return pn[v];
        have[v] = 1;
        double sum = 0, w = 0;
        for (auto k = nptr_[v]; k < nptr_[v + 1]; ++k)
        {
            const Id c = ncell_[k];
            sum += geo_.cell_volume[c] * pe_[c];
            w += geo_.cell_volume[c];
        }
        return pn[v] = w > 0 ? sum / w : 0.0;
    };
    for (Id e = 0; e < cx.ent[1].n_owned; ++e)
    {
        double ne = 0, wsum = 0;
        for (auto k = geo_.ec_ptr[e]; k < geo_.ec_ptr[e + 1]; ++k)
        {
            const Id c = geo_.ec_cell[k];
            if (cx.cell_material[c] == kSolid) continue;
            ne += geo_.ec_weight[k] * ne_[c];
            wsum += geo_.ec_weight[k];
        }
        if (wsum <= 0) continue;
        ne /= wsum;
        Id a, b;
        cx.edge_nodes(e, a, b);
        if (nptr_[a + 1] == nptr_[a] || nptr_[b + 1] == nptr_[b]) continue;
        E[e] -= prm_.pressure_scale * (node_p(b) - node_p(a)) / std::sqrt(ne * ne + prm_.ne_floor * prm_.ne_floor);
    }
}

void MarsSolver::emf_(std::vector<double> &E)
{
    const CellComplex &cx = lm_.mesh.cx;
    E.assign(cx.n_edges(), 0.0);
    emf_asm_.begin(E.data());
    IdealInputs in;
    in.u_cell = uplus_.data();
    in.B_cell_tot = Btot_.data();
    in.B_cell_ind = Bind_.data();
    in.phi_tot = PhiTot_.data();
    in.phi_ind = Phi_.data();
    add_ideal_emf(lm_, geo_, in, cfg_.ideal, E.data());
    if (cfg_.ambipolar) ambipolar_emf_(E.data());
    const Vec3 E_inf = -1.0 * cross(cfg_.u_inf, cfg_.B_imf);
    for (Id e = 0; e < cx.ent[1].n_owned; ++e)
    {
        if (bc_.edge_mode[e] == 1) E[e] = 0.0;                                 // stationary wall
        else if (bc_.edge_mode[e] == 2) E[e] = dot(E_inf, geo_.edge_vec[e]);   // imposed solar-wind electric field
    }
    emf_asm_.finish(E.data());
}

// ------------------------------------------------------------------------------------------------ chemistry
// Cell-local coupled backward-Euler photochemistry (MarsZ5 ApplyPointImplicitChemistry): densities from the network, momenta divided by
// 1 + dt loss/n, energy with the positive/negative parts of the source split (production at 185 K, loss at the local temperature).
void MarsSolver::chemistry_(double dt)
{
    const CellComplex &cx = lm_.mesh.cx;
    const double g = cfg_.sc.gamma;
    const double dt_s = dt * prm_.time_ref;
    const double energy_scale = 1.0e6 * prm_.L * prm_.kB / (prm_.rho_ref * prm_.U * prm_.U * prm_.U * (g - 1.0));
    const double Efloor = prm_.p_floor / (g - 1.0);
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const std::size_t i = c;
        MARS::ChemistryInput input;
        std::array<double, MARS::NumSpecies> floor_cm3{}, estar{}, speed2{};
        double ne_m3 = 0, pe_Pa = 0;
        for (std::size_t n = 0; n < MARS::NumNeutral; ++n) input.neutral_cm3[n] = std::max(neutral_[i * MARS::NumNeutral + n], 0.0);
        for (std::size_t n = 0; n < MARS::NumPhoto; ++n) input.photo_cm3_s[n] = std::max(photo_[i * MARS::NumPhoto + n], 0.0);
        for (int s = 0; s < NS; ++s)
        {
            const double *u = &U_[s][5 * i];
            const double rho = std::max(u[0], prm_.rho_floor);
            const double sp2 = (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / (rho * rho);
            const double kin = 0.5 * rho * sp2;
            const double Et = std::max(u[4], kin + Efloor);
            const double p = std::max((g - 1.0) * (Et - kin), prm_.p_floor);
            const double n_m3 = rho * prm_.rho_ref / MARS::ParticleMass[s];
            input.ion_cm3[s] = n_m3 * 1e-6;
            floor_cm3[s] = prm_.rho_floor * prm_.rho_ref / MARS::ParticleMass[s] * 1e-6;
            input.ion_temperature_K[s] = p * prm_.p_ref / (std::max(n_m3, 1.0) * prm_.kB);
            estar[s] = Et;
            speed2[s] = sp2;
            ne_m3 += n_m3;
            pe_Pa += p * prm_.p_ref;
        }
        input.electron_temperature_K = pe_Pa / (std::max(ne_m3, 1.0) * prm_.kB);

        MARS::PointImplicitChemistryResult chem = MARS::SolvePointImplicitChemistry(input, dt_s, floor_cm3, 16, 1.0e-8);
        for (int sub = 2; !chem.converged && sub <= 16; sub *= 2)
        {
            MARS::ChemistryInput si = input;
            bool ok = true;
            int iters = 0;
            MARS::PointImplicitChemistryResult last;
            MARS::ChemistryRates avg;
            for (int k = 0; k < sub; ++k)
            {
                si.ion_cm3 = k == 0 ? input.ion_cm3 : last.ion_cm3;
                last = MARS::SolvePointImplicitChemistry(si, dt_s / sub, floor_cm3, 16, 1.0e-8);
                iters += last.iterations;
                if (!last.converged) { ok = false; break; }
                for (int s = 0; s < NS; ++s)
                {
                    avg.production_cm3_s[s] += last.rates.production_cm3_s[s] / sub;
                    avg.loss_cm3_s[s] += last.rates.loss_cm3_s[s] / sub;
                }
                avg.electron_production_cm3_s += last.rates.electron_production_cm3_s / sub;
                avg.electron_loss_cm3_s += last.rates.electron_loss_cm3_s / sub;
            }
            if (ok) { chem = last; chem.iterations = iters; chem.rates = avg; }
        }
        if (!chem.converged)
        {
            std::ostringstream m;
            m << "MarsSolver: point-implicit chemistry failed at cell gid " << cx.ent[3].gid[c] << " dt=" << dt;
            throw std::runtime_error(m.str());
        }
        double ne_new = 0;
        for (double v : chem.ion_cm3) ne_new += v;
        for (int s = 0; s < NS; ++s)
        {
            double *u = &U_[s][5 * i];
            const double n_new = chem.ion_cm3[s];
            const double loss_coef = chem.rates.loss_cm3_s[s] / std::max(n_new, 1.0e-300);
            const double mom_den = 1.0 + dt_s * loss_coef;
            for (int d = 1; d < 4; ++d) u[d] /= mom_den;
            const double loss_mass = chem.rates.loss_cm3_s[s] * 1.0e6 * MARS::ParticleMass[s] / prm_.rho_ref * prm_.time_ref;
            const double efrac = n_new / std::max(ne_new, 1.0e-300);
            const double pos = energy_scale * (chem.rates.production_cm3_s[s] * 185.0 + efrac * chem.rates.electron_production_cm3_s * 185.0);
            const double neg = energy_scale * (chem.rates.loss_cm3_s[s] * std::max(input.ion_temperature_K[s], 50.0) +
                                               efrac * chem.rates.electron_loss_cm3_s * std::max(input.electron_temperature_K, 50.0)) +
                               0.5 * loss_mass * speed2[s];
            const double e_den = 1.0 + dt * neg / std::max(estar[s], Efloor);
            u[0] = n_new * 1.0e6 * MARS::ParticleMass[s] / prm_.rho_ref;
            u[4] = (estar[s] + dt * pos) / e_den;
            if (!std::isfinite(u[0]) || !std::isfinite(u[1]) || !std::isfinite(u[2]) || !std::isfinite(u[3]) || !std::isfinite(u[4]) || !(mom_den > 0) || !(e_den > 0))
            {
                std::ostringstream m;
                m << "MarsSolver: invalid chemistry state, species " << name_of(s) << " cell gid " << cx.ent[3].gid[c] << " rho=" << u[0] << " E=" << u[4];
                throw std::runtime_error(m.str());
            }
        }
    }
}

void MarsSolver::floors_()
{
    const CellComplex &cx = lm_.mesh.cx;
    const double g = cfg_.sc.gamma;
    for (int s = 0; s < NS; ++s)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            double *u = &U_[s][5 * S(c)];
            double rho = u[0];
            if (!std::isfinite(rho) || rho < prm_.rho_floor)
            {
                const double scale = prm_.rho_floor / std::max(rho, 1e-30);
                for (int m = 1; m < 4; ++m) u[m] *= scale;
                u[0] = rho = prm_.rho_floor;
            }
            const double ke = 0.5 * (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / rho;
            u[4] = std::max(u[4], ke + prm_.p_floor / (g - 1.0));
        }
}

// ------------------------------------------------------------------------------------------------ one step
double MarsSolver::step(double dt_max)
{
    if (!initialized_) throw std::logic_error("MarsSolver::step before initialize()");
    const CellComplex &cx = lm_.mesh.cx;
    prepare();
    double dt = local_dt_();
    ex_.allreduce(&dt, 1, MPI_MIN);
    dt = std::min(dt, dt_max);
    if (!std::isfinite(dt) || dt <= 0) throw std::runtime_error("MarsSolver: invalid time step");
    fluid_and_sources_(dt);
    std::vector<double> E;
    emf_(E);
    for (int s = 0; s < NS; ++s)
        for (Id c = 0; c < cx.ent[3].n_owned; ++c)
        {
            if (cx.cell_material[c] == kSolid) continue;
            for (int m = 0; m < 5; ++m) U_[s][5 * S(c) + m] += dt * rhs_[s][5 * S(c) + m];
        }
    if (ind_.active()) // Hall term: implicit in the new flux, E is the unknown; wall and inflow edges keep their prescribed values
    {
        std::vector<double> ne(cx.n_cells(), 0.0);
        for (Id c = 0; c < cx.n_cells(); ++c)
            if (cx.cell_material[c] == kFluid) ne[c] = ne_eff_[c];
        ind_.set_frozen_state(Btot_.data(), ne.data());
        std::vector<char> dir(cx.n_edges(), 0);
        for (Id e = 0; e < cx.n_edges(); ++e) dir[e] = bc_.edge_mode[e] != 0;
        std::vector<double> Eimp = E;
        double res = 0;
        gmres_total_ += ind_.solve(Phi_.data(), E.data(), dir, dt, Eimp.data(), &res);
        E.swap(Eimp);
    }
    faraday_update(cx, active_, E.data(), Phi_.data(), dt);
    implicit_sources_(dt);
    chemistry_(dt);
    floors_();
    time_ += dt;
    ++nstep_;
    return dt;
}

// ------------------------------------------------------------------------------------------------ diagnostics, output
MarsSolver::Diag MarsSolver::diagnostics()
{
    const CellComplex &cx = lm_.mesh.cx;
    prepare();
    Diag d;
    d.time = time_; d.step = nstep_;
    double mass[NS], rmin[NS], rmax[NS], pmin[NS], umax[NS], misc[3] = {mhd::div_b(cx, Phi_.data(), PhiAdd_.data()).max_rel, 0, 0};
    for (int s = 0; s < NS; ++s) { mass[s] = 0; rmin[s] = pmin[s] = 1e300; rmax[s] = umax[s] = -1e300; }
    for (Id c = 0; c < cx.ent[3].n_owned; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        for (int s = 0; s < NS; ++s)
        {
            const double *w = &W_[s][5 * S(c)];
            mass[s] += U_[s][5 * S(c)] * geo_.cell_volume[c];
            rmin[s] = std::min(rmin[s], w[0]); rmax[s] = std::max(rmax[s], w[0]);
            pmin[s] = std::min(pmin[s], w[4]);
            umax[s] = std::max(umax[s], std::sqrt(w[1] * w[1] + w[2] * w[2] + w[3] * w[3]));
        }
        misc[1] = std::max(misc[1], norm(Vec3(Btot_[3 * S(c)], Btot_[3 * S(c) + 1], Btot_[3 * S(c) + 2])));
        misc[2] = std::max(misc[2], norm(Vec3(J_[3 * S(c)], J_[3 * S(c) + 1], J_[3 * S(c) + 2])));
    }
    double gm[NS], gmin[2 * NS], gmax[2 * NS], gmisc[3], lo[2 * NS], hi[2 * NS];
    for (int s = 0; s < NS; ++s) { lo[s] = rmin[s]; lo[NS + s] = pmin[s]; hi[s] = rmax[s]; hi[NS + s] = umax[s]; }
    MPI_Allreduce(mass, gm, NS, MPI_DOUBLE, MPI_SUM, comm_);
    MPI_Allreduce(lo, gmin, 2 * NS, MPI_DOUBLE, MPI_MIN, comm_);
    MPI_Allreduce(hi, gmax, 2 * NS, MPI_DOUBLE, MPI_MAX, comm_);
    MPI_Allreduce(misc, gmisc, 3, MPI_DOUBLE, MPI_MAX, comm_);
    for (int s = 0; s < NS; ++s) { d.mass[s] = gm[s]; d.rho_min[s] = gmin[s]; d.p_min[s] = gmin[NS + s]; d.rho_max[s] = gmax[s]; d.u_max[s] = gmax[NS + s]; }
    d.div_b = gmisc[0]; d.b_max = gmisc[1]; d.j_max = gmisc[2];
    // where the time step is set (rank with the smallest limit reports it)
    local_dt_();
    struct { double v; int r; } in{lim_.dt, rank_}, out;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm_);
    double buf[4] = {lim_.x.x, lim_.x.y, lim_.x.z, double(lim_.kind)};
    MPI_Bcast(buf, 4, MPI_DOUBLE, out.r, comm_);
    d.dt_x[0] = buf[0]; d.dt_x[1] = buf[1]; d.dt_x[2] = buf[2]; d.dt_kind = static_cast<int>(buf[3]);
    return d;
}

void MarsSolver::log_(const Diag &d) const
{
    if (rank_ != 0) return;
    std::printf("step %7ld  t=%.6e  dt=%.3e  divB=%.2e  |B|max=%.3e |J|max=%.3e  gmres/step=%.1f\n", d.step, d.time, d.dt, d.div_b, d.b_max, d.j_max,
                nstep_ > 0 ? double(gmres_total_) / double(nstep_) : 0.0);
    std::printf("      dt set by %s at (%.3f %.3f %.3f) r=%.4f\n", d.dt_kind >= NS ? "fast mode" : name_of(std::max(d.dt_kind, 0)), d.dt_x[0], d.dt_x[1], d.dt_x[2],
                std::sqrt(d.dt_x[0] * d.dt_x[0] + d.dt_x[1] * d.dt_x[1] + d.dt_x[2] * d.dt_x[2]));
    for (int s = 0; s < NS; ++s)
        std::printf("      %-4s mass=%.5e  rho=[%.3e, %.3e]  pmin=%.3e  |u|max=%.3e\n", name_of(s), d.mass[s], d.rho_min[s], d.rho_max[s], d.p_min[s], d.u_max[s]);
    std::fflush(stdout);
}

void MarsSolver::run()
{
    if (!initialized_) initialize();
    make_dir(cfg_.output_dir);
    log_(diagnostics());
    if (cfg_.write_vtk) write_vtk();
    double dt = 0;
    while (time_ < cfg_.t_end * (1.0 - 1e-12) && nstep_ < cfg_.max_steps)
    {
        dt = step(cfg_.t_end - time_);
        if (cfg_.log_every > 0 && nstep_ % cfg_.log_every == 0) { Diag d = diagnostics(); d.dt = dt; log_(d); }
        if (cfg_.output_every > 0 && nstep_ % cfg_.output_every == 0 && cfg_.write_vtk) write_vtk();
        if (cfg_.restart_every > 0 && nstep_ % cfg_.restart_every == 0) write_global_fields(cfg_.output_dir + "/fields_final.unsb");
    }
    Diag d = diagnostics();
    d.dt = dt;
    log_(d);
    if (cfg_.write_vtk) write_vtk();
    if (cfg_.write_fields) write_global_fields(cfg_.output_dir + "/fields_final.unsb");
    if (rank_ == 0) std::printf("MarsUns finished normally (%ld steps, t = %.6e)\n", nstep_, time_);
}

void MarsSolver::write_vtk()
{
    prepare();
    const CellComplex &cx = lm_.mesh.cx;
    const Id n = cx.ent[3].n_owned;
    std::vector<io::VtkField> f;
    for (int s = 0; s < NS; ++s) f.push_back({std::string("rho_") + name_of(s), 1, {}});
    for (int s = 0; s < NS; ++s) f.push_back({std::string("p_") + name_of(s), 1, {}});
    for (int s = 0; s < NS; ++s) f.push_back({std::string("u_") + name_of(s), 3, {}});
    f.push_back({"B_induced", 3, {}}); f.push_back({"B_total", 3, {}}); f.push_back({"J", 3, {}}); f.push_back({"u_plus", 3, {}});
    f.push_back({"n_e", 1, {}}); f.push_back({"neutral_CO2", 1, {}}); f.push_back({"neutral_O", 1, {}}); f.push_back({"material", 1, {}});
    for (Id c = 0; c < n; ++c)
    {
        const std::size_t i = c;
        for (int s = 0; s < NS; ++s)
        {
            f[s].values.push_back(W_[s][5 * i]);
            f[NS + s].values.push_back(W_[s][5 * i + 4]);
            for (int a = 0; a < 3; ++a) f[2 * NS + s].values.push_back(W_[s][5 * i + 1 + a]);
        }
        const int o = 3 * NS;
        for (int a = 0; a < 3; ++a)
        {
            f[o].values.push_back(Bind_[3 * i + a]); f[o + 1].values.push_back(Btot_[3 * i + a]);
            f[o + 2].values.push_back(J_[3 * i + a]); f[o + 3].values.push_back(uplus_[3 * i + a]);
        }
        f[o + 4].values.push_back(ne_[i]);
        f[o + 5].values.push_back(neutral_[i * MARS::NumNeutral + MARS::CO2]);
        f[o + 6].values.push_back(neutral_[i * MARS::NumNeutral + MARS::O]);
        f[o + 7].values.push_back(cx.cell_material[c]);
    }
    char base[96], piece[256];
    std::snprintf(base, sizeof base, "sol_%07ld", nstep_);
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

void MarsSolver::write_global_fields(const std::string &path)
{
    prepare();
    const CellComplex &cx = lm_.mesh.cx;
    const std::size_t nc = cx.n_cells();
    std::vector<double> rho(nc, 0.0), p(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c)
        for (int s = 0; s < NS; ++s) { rho[c] += W_[s][5 * c]; p[c] += W_[s][5 * c + 4]; }
    struct Item { std::string name; Dim dim; int ncomp; const double *v; };
    std::vector<Item> items;
    for (int s = 0; s < NS; ++s) items.push_back({std::string("U_") + name_of(s), Dim::Cell, 5, U_[s].data()});
    items.push_back({"B", Dim::Face, 1, Phi_.data()});
    items.push_back({"Badd", Dim::Face, 1, PhiAdd_.data()});
    items.push_back({"rho", Dim::Cell, 1, rho.data()});
    items.push_back({"p", Dim::Cell, 1, p.data()});
    items.push_back({"B_cell", Dim::Cell, 3, Btot_.data()});
    items.push_back({"J_cell", Dim::Cell, 3, J_.data()});
    std::vector<io::FieldRecord> out;
    for (const Item &it : items)
    {
        auto g = gather_entities(comm_, cx.ent[di(it.dim)], it.v, it.ncomp);
        if (rank_ == 0) out.push_back({it.name, it.dim, it.ncomp, std::move(g)});
    }
    if (rank_ == 0)
    {
        out.push_back({"_time", Dim::Node, 1, {time_}});
        io::write_fields(path, 0, out);
        std::printf("wrote %s (t = %.6e)\n", path.c_str(), time_);
    }
}
} // namespace mars
