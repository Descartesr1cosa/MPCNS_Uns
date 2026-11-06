#pragma once
// Initial state and static background field.
#include <vector>

#include "config/Config.h"
#include "uns/geometry/Geometry.h"
#include "uns/parallel/LocalMesh.h"

namespace mercury
{
using namespace mhd;
// Free stream everywhere in the fluid, solid state in solid cells; B_imf as uniform induced flux.
// If cfg.init_fields is set, U_H, U_Na (cells, 5) and B (faces) are taken from that global file by global id.
void set_initial_state(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg, std::vector<double> &UH,
                       std::vector<double> &UNa, std::vector<double> &Phi, double *time_out = nullptr);

// Static dipole: face flux d1(A.dl) of the vector potential (exactly divergence free) and its cell vectors.
void set_background(const uns::LocalMesh &lm, const uns::Geometry &geo, const Config &cfg, std::vector<double> &PhiAdd,
                    std::vector<double> &BaddCell);
} // namespace mercury
