#include "physics/Exosphere.h"

#include <algorithm>
#include <cmath>

namespace mercury
{
using namespace mhd;
using uns::Vec3;

double neutral_sodium_cm3(const Vec3 &x, double Rm)
{
    const double r = uns::norm(x);
    if (!(r > 1.0) || !std::isfinite(r)) return 0.0;
    const double deg = std::atan(1.0) / 45.0;
    const double n_miv0 = 7.84e6, H_miv = 431.0e3;
    const double n_td0 = 8.86e9, H_td = 100.0e3, Ht_td = 15 * deg;
    const double n_psd0 = 4.06e10, H_psd = 232.0e3, Ht_psd = 20 * deg;
    const double n_sp0 = 5.67e6, H_sp = 748.0e3, Ht_sp_day = 15 * deg, Ht_sp_night = 10 * deg;
    const double pi = 4 * std::atan(1.0);
    const Vec3 rh = (1.0 / r) * x;
    const double dr = Rm * (r - 1.0);
    auto bell = [&](double lat, double lon, double Ht) {
        const Vec3 pk(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
        const double a = std::acos(std::max(-1.0, std::min(1.0, dot(rh, pk)))) / Ht;
        return std::exp(-a * a);
    };
    const double n_miv = n_miv0 * std::exp(-dr / H_miv);
    const double n_td = n_td0 * std::exp(-dr / H_td) * bell(0.0, 0.0, Ht_td);
    const double n_psd = n_psd0 * std::exp(-dr / H_psd) * (bell(50 * deg, 0.0, Ht_psd) + bell(-50 * deg, 0.0, Ht_psd));
    const double rad_sp = n_sp0 * std::exp(-dr / H_sp);
    const double n_sp = rad_sp * (bell(80 * deg, 0.0, Ht_sp_day) + bell(-80 * deg, 0.0, Ht_sp_day) + bell(15 * deg, pi, Ht_sp_night) +
                                  bell(-15 * deg, pi, Ht_sp_night));
    double n = (n_miv + n_td + n_psd + n_sp) / (r * r); // [m^-3]
    if (!std::isfinite(n) || n < 0) n = 0;
    return n * 1.0e-6;
}

double photo_rate(const Vec3 &x, double Rm, const ChemistryRates &c)
{
    const double q = (in_shadow(x) ? c.sk2 : c.sk1) * neutral_sodium_cm3(x, Rm);
    return std::isfinite(q) && q > 0 ? q : 0.0;
}
} // namespace mercury
