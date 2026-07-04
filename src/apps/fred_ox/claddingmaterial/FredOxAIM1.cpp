#include "FredOxAIM1.hpp"
#include <cmath>
#include <algorithm>

// Correlations from legacy FRED.f90: clamb, ccp, ctexp, celmod, cpoir, cmhard, ccreep.

namespace fred {

static constexpr double T_REF_AIM1 = 293.15;

double FredOxAIM1::thermalConductivity(double T) const {
    double Tc = T - 273.15;
    return 13.95 + 1.163e-2 * Tc;
}

double FredOxAIM1::heatCapacity(double T) const {
    return 431.0 + 0.177 * T + 8.72e-5 * T * T;
}

double FredOxAIM1::thermalExpansionStrain(double T) const {
    auto eps = [](double t) {
        return -0.2177e-2 + 6.735e-6 * t + 5.12e-9 * t * t
               - 2.248e-12 * t * t * t + 3.933e-16 * t * t * t * t;
    };
    return eps(T) - eps(T_REF_AIM1);
}

double FredOxAIM1::youngsModulus(double T) const {
    double Tc = T - 273.15;
    return (2.027e11 - 8.167e7 * Tc) / 1.0e6;
}

double FredOxAIM1::poissonRatio() const { return 0.289; }

double FredOxAIM1::meyerHardness(double T) const {
    double H = (T <= 893.9203) ? 5.961e9 * std::pow(T, -0.206)
                               : 2.750e28 * std::pow(T, -6.530);
    return std::max(H, 1.0e5);
}

double FredOxAIM1::referenceDensity() const { return m_rho0; }

double FredOxAIM1::creepRate(double T, double sigma) const {
    // Luzzi et al., RdS/PAR2013/022; fixed E=2 MeV, nflux=1e18 n/(cm2·s)
    const double tt = 1.986 * T; // R*T [cal/mol]
    double ccreep = 2.3e14 * std::exp(-8.46e4 / tt) * std::sinh(39.72 * sigma / tt);
    ccreep += 3.2e-24 * 2.0 * 1.0e18 * sigma; // irradiation creep (%/h)
    return ccreep / 100.0 / 3600.0;            // %/h → 1/s
}

} // namespace fred
