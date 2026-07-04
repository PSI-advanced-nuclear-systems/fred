#include "UO2.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

static constexpr double UO2_THEO_DENSITY  = 10960.0; // kg/m3
static constexpr double UO2_MELTING_TEMP  = 3120.0;  // K (pure UO2)
static constexpr double T_REF             = 293.15;  // K

UO2::UO2(double rho0) : m_rho0(rho0) {}

// Philipponneau (1992) at Pu=0, zero burnup, as-fabricated density.
double UO2::thermalConductivity(double T) const {
    // ac at zero burnup and x=0 (stoichiometric)
    double ac  = 1.320 * std::sqrt(0.0093) - 0.091;
    double por = 1.0 - m_rho0 / UO2_THEO_DENSITY;
    return (1.0 / (ac + 2.493e-4 * T) + 88.4e-12 * T * T * T)
           * (1.0 - por) / (1.0 + 2.0 * por);
}

// Popov et al. ORNL/TM-2000/351, pure UO2 tabular data, linear interpolation.
double UO2::heatCapacity(double T) const {
    static const double tTab[] = {
        300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
        1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,
        2500,2600,2700,2800,2900,3000,3100
    };
    static const double cpTab[] = {
        235.51,265.79,282.14,292.21,299.11,304.24,308.32,311.74,
        314.76,317.59,320.44,323.60,327.41,332.31,338.77,347.30,
        358.40,372.54,390.12,411.46,436.78,466.21,499.80,537.50,
        579.17,624.63,673.63,725.88,781.08
    };
    constexpr int N = 29;
    if (T <= tTab[0])   return cpTab[0];
    if (T >= tTab[N-1]) return cpTab[N-1];
    int i = 0;
    while (i < N-1 && T > tTab[i+1]) ++i;
    double f = (T - tTab[i]) / (tTab[i+1] - tTab[i]);
    return cpTab[i] + f * (cpTab[i+1] - cpTab[i]);
}

// MATPRO linear thermal expansion strain (relative to T_REF=293.15 K).
double UO2::thermalExpansionStrain(double T) const {
    auto eps = [](double t) {
        return 1.0e-5 * t - 3.0e-3 + 4.0e-2 * std::exp(-5000.0 / t);
    };
    return eps(T) - eps(T_REF);
}

// MATPRO Young's modulus [MPa] with porosity correction and near-melting softening.
double UO2::youngsModulus(double T, double density) const {
    double por = 1.0 - density / UO2_THEO_DENSITY;
    double E = 2.334e11 * (1.0 - 2.752 * por) * (1.0 - 1.0915e-4 * T);
    const double t0 = UO2_MELTING_TEMP - 150.0;
    if (T > t0) {
        double ratio = (T - t0) / (UO2_MELTING_TEMP - t0);
        double fsoft = std::max(0.0, 1.0 - std::pow(ratio, 12.0));
        E *= fsoft;
    }
    return E / 1.0e6;
}

double UO2::poissonRatio()         const { return 0.276; }
double UO2::referenceDensity()     const { return m_rho0; }
double UO2::theoreticalDensity()   const { return UO2_THEO_DENSITY; }
double UO2::meltingTemperature()   const { return UO2_MELTING_TEMP; }

} // namespace fred
