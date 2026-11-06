#pragma once
// Neutral sodium exosphere (Exner et al. 2020, Table 2: MIV + TD + PSD + SP components) and the Na+ photo-ionisation rate,
// as in the legacy Mercury_Initial::Initial_Na_And_PhotoRate_.  Positions are in planet radii, x > 0 is the dayside.
#include "physics/ChemistryRates.h"
#include "uns/core/Types.h"

namespace mercury
{
using namespace mhd;
// neutral Na density [cm^-3]; zero inside r <= 1
double neutral_sodium_cm3(const uns::Vec3 &x, double L_ref_m);
inline bool in_shadow(const uns::Vec3 &x) { return x.x < 0.0 && x.y * x.y + x.z * x.z < 1.0; }
// production rate q = nu * n_Na  [cm^-3 s^-1]
double photo_rate(const uns::Vec3 &x, double L_ref_m, const ChemistryRates &c);
} // namespace mercury
