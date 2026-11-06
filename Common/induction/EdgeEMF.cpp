#include "induction/EdgeEMF.h"

#include <algorithm>

namespace mhd
{
void EdgeEMF::begin(double *E) const
{
    std::fill(E, E + lm_.mesh.cx.n_edges(), 0.0);
}

double EdgeEMF::weight(uns::Id e, uns::Id c) const
{
    const auto lo = g_.ec_cell.begin() + g_.ec_ptr[e], hi = g_.ec_cell.begin() + g_.ec_ptr[e + 1];
    const auto it = std::lower_bound(lo, hi, c);
    return it != hi && *it == c ? g_.ec_weight[it - g_.ec_cell.begin()] : 0.0;
}
} // namespace mhd
