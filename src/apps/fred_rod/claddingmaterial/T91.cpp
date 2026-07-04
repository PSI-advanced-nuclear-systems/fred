#include "T91.hpp"
#include <cmath>
#include <algorithm>

// All correlations ported from legacy FRED.f90 cmat='t91' branches.
// References are documented in T91.hpp.

namespace fred {

static constexpr double T_REF = 293.15; // K

T91::T91(double rho0) : m_rho0(rho0) {}

// FRED.f90 clamb for cmat='t91': k = 23.71 + 0.01718*Tc - 1.45e-5*Tc^2
// At 300 K (27°C):  k = 23.71 + 0.470 - 0.011 = 24.17 W/(m·K)
// At 823 K (550°C): k = 23.71 + 9.45  - 4.38  = 28.78 W/(m·K)
double T91::thermalConductivity(double T) const {
    double Tc = T - 273.15;
    return 23.71 + 0.01718 * Tc - 1.45e-5 * Tc * Tc;
}

// FRED.f90 ccp for cmat='t91' (identical to AIM1 in legacy code).
double T91::heatCapacity(double T) const {
    return 431.0 + 0.177 * T + 8.72e-5 * T * T;
}

// FRED.f90 ctexp for cmat='t91' (identical polynomial to AIM1).
// Absolute strain: eps(T) = -0.2177e-2 + 6.735e-6*T + 5.12e-9*T^2 - 2.248e-12*T^3 + 3.933e-16*T^4
// Returned = eps(T) - eps(T_REF) so strain is zero at 293.15 K.
double T91::thermalExpansionStrain(double T) const {
    auto eps = [](double t) {
        return -0.2177e-2 + 6.735e-6 * t + 5.12e-9 * t * t
               - 2.248e-12 * t * t * t + 3.933e-16 * t * t * t * t;
    };
    return eps(T) - eps(T_REF);
}

// FRED.f90 celmod for cmat='t91'.
// Two-piece linear in Tc (°C):
//   Tc <= 500°C : E = 2.073e11 - 6.458e7*Tc  Pa
//   Tc >  500°C : E = 2.95e11  - 2.40e8*Tc   Pa
// At Tc=0   (273 K): E = 207.3 GPa
// At Tc=500 (773 K): E = 207.3 - 32.3 = 175.0 GPa  (low-T branch)
//                      = 295.0 - 120.0 = 175.0 GPa  (high-T branch, continuous)
// At Tc=600 (873 K): E = 295.0 - 144.0 = 151.0 GPa
double T91::youngsModulus(double T) const {
    double Tc = T - 273.15;
    double E_Pa;
    if (Tc <= 500.0)
        E_Pa = 2.073e11 - 6.458e7 * Tc;
    else
        E_Pa = 2.95e11  - 2.40e8  * Tc;
    return E_Pa / 1.0e6; // Pa → MPa
}

// FRED.f90 cpoir for cmat='t91'.
double T91::poissonRatio() const { return 0.28; }

// SAS4A correlation for SS Meyer hardness — same formula used for AIM1 and T91
// in FRED.f90 gaphtc (the hardness branch for closed-gap conductance).
// H_M(T<=893.9 K) = 5.961e9 * T^(-0.206)  Pa
// H_M(T> 893.9 K) = 2.750e28 * T^(-6.530)  Pa
// Lower bound: 1e5 Pa to prevent divide-by-zero in gap model.
double T91::meyerHardness(double T) const {
    double H = (T <= 893.9203) ? 5.961e9 * std::pow(T, -0.206)
                                : 2.750e28 * std::pow(T, -6.530);
    return std::max(H, 1.0e5);
}

double T91::referenceDensity() const { return m_rho0; }

} // namespace fred
