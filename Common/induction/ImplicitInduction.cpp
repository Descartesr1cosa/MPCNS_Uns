#include "induction/ImplicitInduction.h"

#include <algorithm>
#include <cmath>

namespace mhd
{
using namespace uns;

ImplicitInduction::ImplicitInduction(const LocalMesh &lm, const Geometry &geo, comm::Exchanger &ex, const Config &cfg)
    : lm_(lm), geo_(geo), ex_(ex), resistive_(cfg.resistive && cfg.sc.inver_Rem > 0), hall_(cfg.hall && cfg.sc.hall_coef > 0), ideal_(cfg.implicit_ideal), ideal_prm_(cfg.ideal), rtol_(cfg.ind_rtol),
      max_it_(cfg.ind_max_it), restart_(std::max(2, cfg.ind_restart)), use_ilu_(cfg.ind_ilu), shift_(cfg.ind_ilu_shift)
{
    const CellComplex &cx = lm.mesh.cx;
    const Id ne = cx.n_edges(), nc = cx.n_cells();
    m_.assign(ne, 0.0);
    if (resistive_)
        for (Id e = 0; e < ne; ++e)
        {
            const double r = norm(geo.edge_mid[e]);
            const double shape = 0.5 * (std::tanh((r - cfg.eta_r_inner) / cfg.eta_width) - std::tanh((r - cfg.eta_r_outer) / cfg.eta_width));
            if (std::isfinite(shape) && shape > 1e-10) m_[e] = cfg.sc.inver_Rem * std::min(1.0, shape) / geo.star1[e];
        }
    taper_.assign(nc, 0.0);
    alpha_.assign(nc, 0.0);
    B0_.assign(3 * static_cast<std::size_t>(nc), 0.0);
    const double mid = 0.5 * (cfg.hall_r_min + cfg.hall_r_max), width = cfg.hall_r_max - cfg.hall_r_min;
    for (Id c = 0; c < nc; ++c)
    {
        if (cx.cell_material[c] == kSolid) continue;
        const double r = norm(geo.cell_center[c]);
        taper_[c] = width <= 0 ? (r >= cfg.hall_r_min ? 1.0 : 0.0) : 0.5 * (1.0 + std::tanh(2.0 * (r - mid) / width));
    }
    hall_coef_ = cfg.sc.hall_coef;

    // cell reconstruction coefficients: B_c = M^-1 sum_f (s n_f) Phi_f  with M = sum_f A_f n n^T
    recon_.assign(3 * static_cast<std::size_t>(cx.c2f.nnz()), 0.0);
    for (Id c = 0; c < nc; ++c)
    {
        double M[3][3] = {{0}};
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            const double A = geo.face_area_mag[f];
            if (A <= 0) continue;
            const Vec3 n = (1.0 / A) * geo.face_area[f];
            const double nv[3] = {n.x, n.y, n.z};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) M[i][j] += A * nv[i] * nv[j];
        }
        const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        const double sc = std::max({M[0][0], M[1][1], M[2][2], 1e-300});
        if (std::abs(det) < 1e-12 * sc * sc * sc) continue;
        double I[3][3];
        I[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det; I[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / det; I[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
        I[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) / det; I[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det; I[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / det;
        I[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det; I[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / det; I[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
        for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
        {
            const Id f = cx.c2f.idx[k];
            const double A = geo.face_area_mag[f];
            if (A <= 0) continue;
            const Vec3 n = (1.0 / A) * geo.face_area[f]; // outward normal x outward flux = canonical normal x canonical flux
            const double nv[3] = {n.x, n.y, n.z};
            for (int i = 0; i < 3; ++i) recon_[3 * static_cast<std::size_t>(k) + i] = I[i][0] * nv[0] + I[i][1] * nv[1] + I[i][2] * nv[2];
        }
    }

    // owned cell -> (edge, w_{e,c} t_e)
    const Id no = cx.ent[3].n_owned;
    ce_ptr_.assign(no + 1, 0);
    std::vector<Id> es;
    for (Id c = 0; c < no; ++c)
    {
        es.clear();
        if (cx.cell_material[c] == kFluid)
            for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
                for (auto ke = cx.f2e.begin(cx.c2f.idx[kf]); ke < cx.f2e.end(cx.c2f.idx[kf]); ++ke) es.push_back(cx.f2e.idx[ke]);
        std::sort(es.begin(), es.end());
        es.erase(std::unique(es.begin(), es.end()), es.end());
        for (Id e : es)
        {
            const auto lo = geo.ec_cell.begin() + geo.ec_ptr[e], hi = geo.ec_cell.begin() + geo.ec_ptr[e + 1];
            const auto it = std::lower_bound(lo, hi, c);
            if (it == hi || *it != c) continue;
            ce_edge_.push_back(e);
            ce_wt_.push_back(geo.ec_weight[it - geo.ec_cell.begin()] * geo.edge_vec[e]);
        }
        ce_ptr_[c + 1] = static_cast<std::int64_t>(ce_edge_.size());
    }
    d2_.assign(ne, 0.0);
    for (Id f = 0; f < cx.ent[2].n_owned; ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) d2_[cx.f2e.idx[k]] += geo.star2[f];
    ex_.assemble_edges(d2_.data());
    phi_.assign(cx.n_faces(), 0.0);
    bcell_.assign(3 * static_cast<std::size_t>(nc), 0.0);
    jcell_.assign(3 * static_cast<std::size_t>(nc), 0.0);
    part_.assign(3 * static_cast<std::size_t>(ne), 0.0);
    ide_.assign(ne, 0.0);
    tmp_.assign(ne, 0.0);
    if (use_ilu_ && active() && !ideal_) build_precond_pattern();
}

void ImplicitInduction::set_frozen_state(const double *B0, const double *ne_eff)
{
    const CellComplex &cx = lm_.mesh.cx;
    for (std::size_t i = 0; i < B0_.size(); ++i) B0_[i] = B0[i];
    for (Id c = 0; c < cx.n_cells(); ++c)
        alpha_[c] = (hall_ && ne_eff[c] > 0) ? taper_[c] * hall_coef_ / ne_eff[c] : 0.0;
}

void ImplicitInduction::set_ideal_state(const double *u_cell, const std::vector<char> &mask)
{
    u_ideal_ = u_cell;
    ideal_mask_ = mask;
}

void ImplicitInduction::ideal_linear(const double *phi, double *out)
{
    const bool h = hall_, r = resistive_;
    hall_ = false; resistive_ = false; // reuse apply_F with only the ideal part switched on
    apply_F(phi, out);
    hall_ = h; resistive_ = r;
}

void ImplicitInduction::apply_F(const double *phi, double *out)
{
    const CellComplex &cx = lm_.mesh.cx;
    const Id ne = cx.n_edges(), nc = cx.n_cells(), no = cx.ent[3].n_owned;
    std::fill(part_.begin(), part_.end(), 0.0);
    // current flux through the dual faces: d1^T S2 Phi, each face once (by its owner)
    if (resistive_)
        for (Id f = 0; f < cx.ent[2].n_owned; ++f)
        {
            const double v = geo_.star2[f] * phi[f];
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) part_[3 * static_cast<std::size_t>(cx.f2e.idx[k])] += cx.f2e.sgn[k] * v;
        }
    if (hall_ || ideal_)
        for (Id c = 0; c < nc; ++c) // cell field from the face fluxes
        {
            Vec3 b;
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
            {
                const double v = phi[cx.c2f.idx[k]];
                b.x += recon_[3 * static_cast<std::size_t>(k)] * v; b.y += recon_[3 * static_cast<std::size_t>(k) + 1] * v; b.z += recon_[3 * static_cast<std::size_t>(k) + 2] * v;
            }
            bcell_[3 * c] = b.x; bcell_[3 * c + 1] = b.y; bcell_[3 * c + 2] = b.z;
        }
    if (ideal_ && u_ideal_)
    {
        IdealInputs in;
        in.u_cell = u_ideal_; in.B_cell_tot = bcell_.data(); in.B_cell_ind = bcell_.data(); in.phi_tot = phi; in.phi_ind = phi;
        std::fill(ide_.begin(), ide_.end(), 0.0);
        add_ideal_emf(lm_, geo_, in, ideal_prm_, ide_.data());
        for (Id e = 0; e < ne; ++e) part_[3 * static_cast<std::size_t>(e) + 2] = ideal_mask_[e] ? ide_[e] : 0.0;
    }
    if (hall_)
    {
        for (Id c = 0; c < no; ++c)
        {
            if (alpha_[c] == 0.0) continue;
            Vec3 sum;
            const Vec3 bc(bcell_[3 * c], bcell_[3 * c + 1], bcell_[3 * c + 2]);
            for (auto k = cx.c2f.begin(c); k < cx.c2f.end(c); ++k)
            {
                const Id f = cx.c2f.idx[k];
                Vec3 bf = bc;
                if (cx.f2c.count(f) == 2)
                {
                    const Id o = cx.f2c.idx[cx.f2c.begin(f)] == c ? cx.f2c.idx[cx.f2c.begin(f) + 1] : cx.f2c.idx[cx.f2c.begin(f)];
                    bf = 0.5 * (bc + Vec3(bcell_[3 * o], bcell_[3 * o + 1], bcell_[3 * o + 2]));
                }
                sum += cross(static_cast<double>(cx.c2f.sgn[k]) * geo_.face_area[f], bf);
            }
            const Vec3 J = (1.0 / geo_.cell_volume[c]) * sum;
            const Vec3 g = alpha_[c] * cross(J, Vec3(B0_[3 * c], B0_[3 * c + 1], B0_[3 * c + 2]));
            for (auto q = ce_ptr_[c]; q < ce_ptr_[c + 1]; ++q) part_[3 * static_cast<std::size_t>(ce_edge_[q]) + 1] += dot(g, ce_wt_[q]);
        }
    }
    ex_.assemble_edges(part_.data(), 3);
    for (Id e = 0; e < ne; ++e) out[e] = m_[e] * part_[3 * static_cast<std::size_t>(e)] + part_[3 * static_cast<std::size_t>(e) + 1] + part_[3 * static_cast<std::size_t>(e) + 2];
}

void ImplicitInduction::apply_A(const std::vector<double> &x, std::vector<double> &y, double dt)
{
    const CellComplex &cx = lm_.mesh.cx;
    for (Id f = 0; f < cx.n_faces(); ++f) // Phi' = d1 x on every local face (all its edges are local)
    {
        double s = 0;
        if (dte_) for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) s += cx.f2e.sgn[k] * dte_[cx.f2e.idx[k]] * x[cx.f2e.idx[k]];
        else for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) s += cx.f2e.sgn[k] * x[cx.f2e.idx[k]];
        phi_[f] = s;
    }
    apply_F(phi_.data(), tmp_.data());
    if (dte_) for (std::size_t e = 0; e < y.size(); ++e) y[e] = x[e] + tmp_[e]; // d1 already carries dt_e
    else for (std::size_t e = 0; e < y.size(); ++e) y[e] = x[e] + dt * tmp_[e];
}

double ImplicitInduction::dot_owned(const std::vector<double> &a, const std::vector<double> &b)
{
    double s = 0;
    const Id no = lm_.mesh.cx.ent[1].n_owned;
    for (Id e = 0; e < no; ++e) s += a[e] * b[e];
    ex_.allreduce(&s, 1, MPI_SUM);
    return s;
}

void ImplicitInduction::build_precond_pattern()
{
    const CellComplex &cx = lm_.mesh.cx;
    const Id no = cx.ent[1].n_owned;
    pc_ptr_.assign(no + 1, 0);
    pc_col_.clear(); pc_coef_.clear();
    std::vector<std::pair<Id, double>> row;
    for (Id e = 0; e < no; ++e)
    {
        row.clear();
        row.emplace_back(e, 0.0); // the diagonal is always present
        for (auto kf = cx.e2f.begin(e); kf < cx.e2f.end(e); ++kf)
        {
            const Id f = cx.e2f.idx[kf];
            const double se = cx.e2f.sgn[kf];
            for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k)
            {
                const Id e2 = cx.f2e.idx[k];
                if (e2 >= no) continue; // other ranks: block Jacobi
                row.emplace_back(e2, m_[e] * se * cx.f2e.sgn[k] * geo_.star2[f]);
            }
        }
        std::sort(row.begin(), row.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
        for (std::size_t i = 0; i < row.size();)
        {
            std::size_t j = i;
            double sum = 0;
            for (; j < row.size() && row[j].first == row[i].first; ++j) sum += row[j].second;
            pc_col_.push_back(row[i].first);
            pc_coef_.push_back(sum);
            i = j;
        }
        pc_ptr_[e + 1] = static_cast<std::int64_t>(pc_col_.size());
    }
    pc_val_.assign(pc_col_.size(), 0.0);
    pc_diag_.assign(no, 0);
    for (Id e = 0; e < no; ++e)
        for (auto k = pc_ptr_[e]; k < pc_ptr_[e + 1]; ++k) if (pc_col_[k] == e) pc_diag_[e] = k;
}

void ImplicitInduction::factor_precond(double dt, const std::vector<char> &dir)
{
    const Id no = lm_.mesh.cx.ent[1].n_owned;
    for (Id e = 0; e < no; ++e)
        for (auto k = pc_ptr_[e]; k < pc_ptr_[e + 1]; ++k)
        {
            const Id c = pc_col_[k];
            // Dirichlet unknowns are removed: identity row and no coupling into them
            pc_val_[k] = (dir[e] || dir[c]) ? (c == e ? 1.0 : 0.0) : (c == e ? 1.0 : 0.0) + dt * pc_coef_[k];
            if (c == e && !dir[e]) pc_val_[k] += shift_ * dt * pc_coef_[k]; // relaxed diagonal: keeps the incomplete factors stable
        }
    std::vector<std::int64_t> pos(no, -1);
    for (Id i = 0; i < no; ++i) // ILU(0), IKJ variant
    {
        for (auto k = pc_ptr_[i]; k < pc_ptr_[i + 1]; ++k) pos[pc_col_[k]] = k;
        for (auto kk = pc_ptr_[i]; kk < pc_ptr_[i + 1] && pc_col_[kk] < i; ++kk)
        {
            const Id k = pc_col_[kk];
            pc_val_[kk] /= pc_val_[pc_diag_[k]];
            const double lik = pc_val_[kk];
            for (auto jj = pc_diag_[k] + 1; jj < pc_ptr_[k + 1]; ++jj)
            {
                const std::int64_t p = pos[pc_col_[jj]];
                if (p >= 0) pc_val_[p] -= lik * pc_val_[jj];
            }
        }
        if (!(std::abs(pc_val_[pc_diag_[i]]) > 1e-12)) pc_val_[pc_diag_[i]] = 1.0; // breakdown guard
        for (auto k = pc_ptr_[i]; k < pc_ptr_[i + 1]; ++k) pos[pc_col_[k]] = -1;
    }
}

void ImplicitInduction::apply_precond(const std::vector<double> &v, std::vector<double> &z, const std::vector<char> &dir)
{
    const Id no = lm_.mesh.cx.ent[1].n_owned, ne = lm_.mesh.cx.n_edges();
    for (Id i = 0; i < no; ++i) // L y = v (unit lower)
    {
        double s = dir[i] ? 0.0 : v[i];
        for (auto k = pc_ptr_[i]; k < pc_diag_[i]; ++k) s -= pc_val_[k] * z[pc_col_[k]];
        z[i] = s;
    }
    for (Id i = no; i-- > 0;) // U z = y
    {
        double s = z[i];
        for (auto k = pc_diag_[i] + 1; k < pc_ptr_[i + 1]; ++k) s -= pc_val_[k] * z[pc_col_[k]];
        z[i] = dir[i] ? 0.0 : s / pc_val_[pc_diag_[i]];
    }
    for (Id e = no; e < ne; ++e) z[e] = 0.0;
    ex_.forward(Dim::Edge, z.data(), 1); // refresh the alias copies from the owners
}

int ImplicitInduction::solve(const double *Phi_n, const double *E_expl, const std::vector<char> &dir, double dt, double *E, double *rel_res, const double *dt_edge)
{
    dte_ = dt_edge;
    const CellComplex &cx = lm_.mesh.cx;
    const Id ne = cx.n_edges();
    using Vecd = std::vector<double>;
    Vecd x0(E_expl, E_expl + ne), b(ne), r(ne), Ax(ne);
    if (use_guess_ && corr_.size() == static_cast<std::size_t>(ne)) // warm start: the implicit correction of the previous step
        for (Id e = 0; e < ne; ++e) if (!dir[e]) x0[e] += corr_[e];
    // b = E_expl + F(Phi^n)
    std::vector<double> Fn(ne);
    apply_F(Phi_n, Fn.data());
    for (Id e = 0; e < ne; ++e) b[e] = E_expl[e] + Fn[e];
    apply_A(x0, Ax, dt);
    for (Id e = 0; e < ne; ++e) r[e] = dir[e] ? 0.0 : b[e] - Ax[e];
    // Jacobi diagonal:  1 + dt m_e d2_e  (the Hall part has no useful diagonal)
    Vecd inv_diag(ne);
    for (Id e = 0; e < ne; ++e) inv_diag[e] = 1.0 / (1.0 + (dte_ ? dte_[e] : dt) * m_[e] * d2_[e]);
    const bool ilu = use_ilu_ && !pc_ptr_.empty() && !dte_; // ILU is built for a scalar dt
    if (ilu) factor_precond(dt, dir);

    const double bnorm = std::sqrt(std::max(dot_owned(b, b), 1e-300));
    double beta = std::sqrt(dot_owned(r, r));
    int total = 0;
    std::vector<double> delta(ne, 0.0);
    const int m = restart_;
    std::vector<Vecd> V(m + 1, Vecd(ne));
    std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
    std::vector<double> cs(m), sn(m), g(m + 1), yv(m);
    Vecd w(ne), z(ne);
    while (beta > rtol_ * bnorm && total < max_it_)
    {
        for (Id e = 0; e < ne; ++e) V[0][e] = r[e] / beta;
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;
        int j = 0;
        for (; j < m && total < max_it_; ++j, ++total)
        {
            if (ilu) apply_precond(V[j], z, dir);
            else for (Id e = 0; e < ne; ++e) z[e] = dir[e] ? 0.0 : V[j][e] * inv_diag[e];
            apply_A(z, w, dt);
            for (Id e = 0; e < ne; ++e) if (dir[e]) w[e] = 0.0;
            for (int i = 0; i <= j; ++i) // modified Gram-Schmidt
            {
                H[i][j] = dot_owned(w, V[i]);
                for (Id e = 0; e < ne; ++e) w[e] -= H[i][j] * V[i][e];
            }
            H[j + 1][j] = std::sqrt(dot_owned(w, w));
            if (H[j + 1][j] > 1e-300) for (Id e = 0; e < ne; ++e) V[j + 1][e] = w[e] / H[j + 1][j];
            for (int i = 0; i < j; ++i)
            {
                const double t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = t;
            }
            const double den = std::hypot(H[j][j], H[j + 1][j]);
            cs[j] = den > 0 ? H[j][j] / den : 1.0; sn[j] = den > 0 ? H[j + 1][j] / den : 0.0;
            H[j][j] = den; H[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j]; g[j] = cs[j] * g[j];
            if (std::abs(g[j + 1]) <= rtol_ * bnorm) { ++j; ++total; break; }
        }
        for (int i = j - 1; i >= 0; --i) // back substitution
        {
            double s = g[i];
            for (int k = i + 1; k < j; ++k) s -= H[i][k] * yv[k];
            yv[i] = s / H[i][i];
        }
        {
            Vecd comb(ne, 0.0);
            for (Id e = 0; e < ne; ++e)
            {
                double s = 0;
                for (int i = 0; i < j; ++i) s += yv[i] * V[i][e];
                comb[e] = s;
            }
            if (ilu) { apply_precond(comb, z, dir); for (Id e = 0; e < ne; ++e) delta[e] += z[e]; }
            else for (Id e = 0; e < ne; ++e) delta[e] += dir[e] ? 0.0 : comb[e] * inv_diag[e];
        }
        // true residual
        Vecd xt(ne);
        for (Id e = 0; e < ne; ++e) xt[e] = x0[e] + delta[e];
        apply_A(xt, Ax, dt);
        for (Id e = 0; e < ne; ++e) r[e] = dir[e] ? 0.0 : b[e] - Ax[e];
        beta = std::sqrt(dot_owned(r, r));
    }
    for (Id e = 0; e < ne; ++e) E[e] = x0[e] + delta[e];
    corr_.assign(ne, 0.0);
    for (Id e = 0; e < ne; ++e) corr_[e] = dir[e] ? 0.0 : E[e] - E_expl[e];
    if (rel_res) *rel_res = beta / bnorm;
    dte_ = nullptr;
    return total;
}
} // namespace mhd
