#pragma once
// Electron-pressure (ambipolar) effects of the legacy model, p_e = p_H + p_Na:
//   force on the ions:      S_mom,s = - w_s grad p_e ,  S_E,s = - w_s u_s . grad p_e         (see fluid sources)
//   electric field:         E = - pressure_coef * grad p_e / n_e,eff                          (an edge 1-cochain)
// The edge form is the exact d0 of a nodal pressure divided by an edge-averaged density, so it adds a curl only where
// grad n_e x grad p_e does not vanish (Biermann term) and goes through the same edge assembly as every other EMF.
#include <vector>

#include "config/Config.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
using namespace mhd;
// grad p_e of the owned fluid cells (3 per cell, zero elsewhere) by Gauss' theorem with face-averaged pressure.
void electron_pressure_gradient(const uns::LocalMesh &lm, const uns::Geometry &geo, const double *WH, const double *WNa, double *dpe);

class AmbipolarEMF
{
public:
    AmbipolarEMF(const uns::LocalMesh &lm, const uns::Geometry &geo);
    // adds the ambipolar EMF to the owned edges of E (that touch fluid cells)
    void add(const Config &cfg, const double *WH, const double *WNa, double *E) const;
private:
    const uns::LocalMesh &lm_;
    const uns::Geometry &geo_;
    std::vector<std::int64_t> nptr_;   // node -> fluid cells
    std::vector<uns::Id> ncell_;
};
} // namespace mercury
