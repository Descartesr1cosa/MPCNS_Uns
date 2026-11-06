#include "induction/Resistive.h"

#include <algorithm>
#include <cmath>

namespace mercury
{
using namespace mhd;
using namespace uns;

ResistiveSolver::ResistiveSolver(const LocalMesh &lm, const Geometry &geo, comm::Exchanger &ex, const Config &cfg)
    : lm_(lm), geo_(geo), ex_(ex), rtol_(cfg.res_rtol), max_it_(cfg.res_max_it)
{
    const CellComplex &cx = lm.mesh.cx;
    const Id ne = cx.n_edges();
    in_r_.assign(ne, 0);
    minv_.assign(ne, 0.0);
    if (!cfg.resistive || cfg.sc.inver_Rem <= 0) return;
    for (Id e = 0; e < ne; ++e)
    {
        const double r = norm(geo.edge_mid[e]);
        const double shape = 0.5 * (std::tanh((r - cfg.eta_r_inner) / cfg.eta_width) - std::tanh((r - cfg.eta_r_outer) / cfg.eta_width));
        if (!std::isfinite(shape) || shape <= 1e-10) continue; // numerically negligible tail
        in_r_[e] = 1;
        minv_[e] = geo.star1[e] / (cfg.sc.inver_Rem * std::min(1.0, shape));
        if (e < cx.ent[1].n_owned) ++n_owned_res_;
    }
    double g = static_cast<double>(n_owned_res_);
    ex_.allreduce(&g, 1, MPI_SUM);
    enabled_ = g > 0;
    for (auto *v : {&x_, &b_, &r_, &z_, &p_, &q_}) v->assign(ne, 0.0);
    d2_.assign(ne, 0.0);
    if (!enabled_) return;
    // diagonal of D^T S2 D, assembled over the owned faces
    for (Id f = 0; f < cx.ent[2].n_owned; ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) d2_[cx.f2e.idx[k]] += geo.star2[f];
    ex_.assemble_edges(d2_.data());
}

double ResistiveSolver::dot(const std::vector<double> &a, const std::vector<double> &b)
{
    double s = 0;
    const Id no = lm_.mesh.cx.ent[1].n_owned;
    for (Id e = 0; e < no; ++e) s += a[e] * b[e];
    ex_.allreduce(&s, 1, MPI_SUM);
    return s;
}

// out_e = sum_f d1[f,e] star2_f (sum_e' d1[f,e'] x_e'), faces counted once (by their owner), then assembled through the edge owners
void ResistiveSolver::dt_s2_d(const double *x, double *out) const
{
    const CellComplex &cx = lm_.mesh.cx;
    std::fill(out, out + cx.n_edges(), 0.0);
    for (Id f = 0; f < cx.ent[2].n_owned; ++f)
    {
        double y = 0;
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) y += cx.f2e.sgn[k] * x[cx.f2e.idx[k]];
        y *= geo_.star2[f];
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) out[cx.f2e.idx[k]] += cx.f2e.sgn[k] * y;
    }
    ex_.assemble_edges(out);
}

void ResistiveSolver::apply(const std::vector<double> &x, std::vector<double> &y, double dt)
{
    dt_s2_d(x.data(), y.data());
    for (std::size_t e = 0; e < y.size(); ++e) y[e] = in_r_[e] ? minv_[e] * x[e] + dt * y[e] : 0.0;
}

int ResistiveSolver::advance(double *Phi, double dt, const std::vector<char> &active)
{
    last_it_ = 0;
    if (!enabled_) return 0;
    const CellComplex &cx = lm_.mesh.cx;
    const Id ne = cx.n_edges();
    // right-hand side  b = D^T S2 Phi  on R   (owned faces)
    std::fill(b_.begin(), b_.end(), 0.0);
    for (Id f = 0; f < cx.ent[2].n_owned; ++f)
    {
        const double v = geo_.star2[f] * Phi[f];
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) b_[cx.f2e.idx[k]] += cx.f2e.sgn[k] * v;
    }
    ex_.assemble_edges(b_.data());
    for (Id e = 0; e < ne; ++e) if (!in_r_[e]) b_[e] = 0.0;

    // Jacobi-preconditioned CG, x0 = 0
    std::fill(x_.begin(), x_.end(), 0.0);
    r_ = b_;
    auto precond = [&](const std::vector<double> &r, std::vector<double> &z) {
        for (Id e = 0; e < ne; ++e) z[e] = in_r_[e] ? r[e] / (minv_[e] + dt * d2_[e]) : 0.0;
    };
    precond(r_, z_);
    p_ = z_;
    double rz = dot(r_, z_);
    const double bnorm = std::sqrt(std::max(dot(b_, b_), 1e-300));
    for (int it = 0; it < max_it_; ++it)
    {
        if (std::sqrt(dot(r_, r_)) <= rtol_ * bnorm) break;
        apply(p_, q_, dt);
        const double alpha = rz / dot(p_, q_);
        for (Id e = 0; e < ne; ++e) { x_[e] += alpha * p_[e]; r_[e] -= alpha * q_[e]; }
        precond(r_, z_);
        const double rz2 = dot(r_, z_);
        const double beta = rz2 / rz;
        rz = rz2;
        for (Id e = 0; e < ne; ++e) p_[e] = z_[e] + beta * p_[e];
        ++last_it_;
    }
    // Phi <- Phi - dt d1 E_res on the faces of owned cells
    for (Id f = 0; f < cx.n_faces(); ++f)
    {
        if (!active[f]) continue;
        double c = 0;
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) c += cx.f2e.sgn[k] * x_[cx.f2e.idx[k]];
        Phi[f] -= dt * c;
    }
    return last_it_;
}
} // namespace mercury
