#include "initial/Initial.h"

#include <stdexcept>

#include "physics/Euler.h"
#include "uns/geometry/Reconstruct.h"
#include "uns/io/MeshIO.h"

namespace mercury
{
using namespace mhd;
using namespace uns;

void set_initial_state(const LocalMesh &lm, const Geometry &geo, const Config &cfg, std::vector<double> &UH, std::vector<double> &UNa,
                       std::vector<double> &Phi, double *time_out)
{
    const CellComplex &cx = lm.mesh.cx;
    const std::size_t nc = cx.n_cells();
    UH.assign(5 * nc, 0.0);
    UNa.assign(5 * nc, 0.0);
    Phi.assign(cx.n_faces(), 0.0);
    auto fill = [&](std::vector<double> &U, const SpeciesState &st, const Vec3 &u) {
        Prim w{st.rho, u.x, u.y, u.z, st.p};
        for (std::size_t c = 0; c < nc; ++c) to_cons(w, cfg.sc.gamma, &U[5 * c]);
    };
    fill(UH, cfg.fsH, cfg.u_inf);
    fill(UNa, cfg.fsNa, cfg.u_Na());
    const Vec3 zero;
    for (std::size_t c = 0; c < nc; ++c)
        if (cx.cell_material[c] == kSolid)
        {
            Prim w{cfg.solidState.rho, 0, 0, 0, cfg.solidState.p};
            to_cons(w, cfg.sc.gamma, &UH[5 * c]);
            to_cons(w, cfg.sc.gamma, &UNa[5 * c]);
        }
    // induced field starts at zero (the IMF is part of the static background, see set_background)

    if (!cfg.init_fields.empty())
    {
        std::uint64_t hash = 0;
        const auto fields = io::read_fields(cfg.init_fields, &hash);
        bool gotH = false, gotNa = false, gotB = false;
        for (const auto &fr : fields)
        {
            if ((fr.name == "U_H" || fr.name == "U_Na") && fr.dim == Dim::Cell && fr.ncomp == 5)
            {
                auto &U = fr.name == "U_H" ? UH : UNa;
                for (std::size_t c = 0; c < nc; ++c)
                    for (int m = 0; m < 5; ++m) U[5 * c + m] = fr.values.at(static_cast<std::size_t>(cx.ent[3].gid[c]) * 5 + m);
                (fr.name == "U_H" ? gotH : gotNa) = true;
            }
            else if (fr.name == "_time" && !fr.values.empty())
            {
                if (time_out) *time_out = fr.values[0];
            }
            else if (fr.name == "B" && fr.dim == Dim::Face)
            {
                for (Id f = 0; f < cx.n_faces(); ++f) Phi[f] = fr.values.at(static_cast<std::size_t>(cx.ent[2].gid[f]));
                gotB = true;
            }
        }
        if (!gotH || !gotNa || !gotB) throw std::runtime_error("init_fields '" + cfg.init_fields + "' lacks U_H, U_Na or B");
    }
}

void set_background(const LocalMesh &lm, const Geometry &geo, const Config &cfg, std::vector<double> &PhiAdd, std::vector<double> &BaddCell)
{
    const CellComplex &cx = lm.mesh.cx;
    PhiAdd.assign(cx.n_faces(), 0.0);
    BaddCell.assign(3 * static_cast<std::size_t>(cx.n_cells()), 0.0);
    for (Id f = 0; f < cx.n_faces(); ++f) PhiAdd[f] = dot(cfg.B_imf, geo.face_area[f]); // uniform IMF
    std::vector<double> A(cx.n_edges(), 0.0);
    for (Id e = 0; cfg.dipole && e < cx.n_edges(); ++e)
    {
        const Vec3 r = geo.edge_mid[e] - cfg.dipole_c;
        const double d2 = dot(r, r) + cfg.dipole_soft * cfg.dipole_soft;
        const Vec3 a = (1.0 / (d2 * std::sqrt(d2))) * cross(cfg.dipole_m, r);
        A[e] = dot(a, geo.edge_vec[e]);
    }
    for (Id f = 0; cfg.dipole && f < cx.n_faces(); ++f)
        for (auto k = cx.f2e.begin(f); k < cx.f2e.end(f); ++k) PhiAdd[f] += cx.f2e.sgn[k] * A[cx.f2e.idx[k]];
    const auto b = reconstruct_cell_vector(cx, geo, PhiAdd.data(), cx.n_cells());
    for (Id c = 0; c < cx.n_cells(); ++c) { BaddCell[3 * c] = b[c].x; BaddCell[3 * c + 1] = b[c].y; BaddCell[3 * c + 2] = b[c].z; }
}
} // namespace mercury
