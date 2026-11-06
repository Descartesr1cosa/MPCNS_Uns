#pragma once
// Rates of the Mercury Na+ source terms carried by the run configuration.
namespace mhd
{
struct ChemistryRates
{
    double sk1 = 5.0e-5;  // photo-ionisation / drag frequency in sunlight  [1/s]
    double sk2 = 1.0e-5;  // ... in the geometric shadow                    [1/s]
    double Tn0 = 1000.0;  // neutral temperature [K]
};
} // namespace mhd
