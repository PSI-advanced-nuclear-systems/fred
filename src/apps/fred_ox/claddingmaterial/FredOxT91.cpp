#include "FredOxT91.hpp"
#include <cmath>
#include <algorithm>

// Correlations from legacy FRED.f90 cmat='t91' branches (clamb, ccp, ctexp, celmod, cpoir).

namespace fred {

static constexpr double T_REF_T91 = 293.15;

double FredOxT91::thermalConductivity(double T) const {
    double Tc = T - 273.15;
    return 23.71 + 0.01718 * Tc - 1.45e-5 * Tc * Tc;
}

double FredOxT91::heatCapacity(double T) const {
    return 431.0 + 0.177 * T + 8.72e-5 * T * T;
}

double FredOxT91::thermalExpansionStrain(double T) const {
    auto eps = [](double t) {
        return -0.2177e-2 + 6.735e-6 * t + 5.12e-9 * t * t
               - 2.248e-12 * t * t * t + 3.933e-16 * t * t * t * t;
    };
    return eps(T) - eps(T_REF_T91);
}

double FredOxT91::youngsModulus(double T) const {
    double Tc = T - 273.15;
    double E_Pa = (Tc <= 500.0) ? (2.073e11 - 6.458e7 * Tc)
                                 : (2.95e11  - 2.40e8  * Tc);
    return E_Pa / 1.0e6;
}

double FredOxT91::poissonRatio() const { return 0.28; }

double FredOxT91::meyerHardness(double T) const {
    double H = (T <= 893.9203) ? 5.961e9 * std::pow(T, -0.206)
                               : 2.750e28 * std::pow(T, -6.530);
    return std::max(H, 1.0e5);
}

double FredOxT91::referenceDensity() const { return m_rho0; }

} // namespace fred
