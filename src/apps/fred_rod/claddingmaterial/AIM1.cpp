#include "AIM1.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

static constexpr double T_REF = 293.15;

AIM1::AIM1(double rho0) : m_rho0(rho0) {}

double AIM1::thermalConductivity(double T) const {
    double Tc = T - 273.15;
    return 13.95 + 1.163e-2 * Tc;
}

double AIM1::heatCapacity(double T) const {
    return 431.0 + 0.177 * T + 8.72e-5 * T * T;
}

double AIM1::thermalExpansionStrain(double T) const {
    auto eps = [](double t) {
        return -0.2177e-2 + 6.735e-6 * t + 5.12e-9 * t * t
               - 2.248e-12 * t * t * t + 3.933e-16 * t * t * t * t;
    };
    return eps(T) - eps(T_REF);
}

double AIM1::youngsModulus(double T) const {
    double Tc = T - 273.15;
    return (2.027e11 - 8.167e7 * Tc) / 1.0e6;
}

double AIM1::poissonRatio() const { return 0.289; }

double AIM1::meyerHardness(double T) const {
    double H = (T <= 893.9203) ? 5.961e9 * std::pow(T, -0.206)
                               : 2.750e28 * std::pow(T, -6.530);
    return std::max(H, 1.0e5);
}

double AIM1::referenceDensity() const { return m_rho0; }

} // namespace fred
