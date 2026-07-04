#include "MOX.hpp"
#include <cmath>
#include <algorithm>

// All correlations ported from legacy FRED.f90 functions: flamb, fcp, ftexp, felmod, fpoir.
// References are documented in MOX.hpp.

namespace fred {

static constexpr double T_REF = 293.15; // K

// Theoretical density formula from FRED.f90 felmod:
//   tden = 11460*Pu + 10960*(1-Pu)  [kg/m3]
// At Pu=0   → 10960 (UO2 TD)
// At Pu=0.15 → 11460*0.15 + 10960*0.85 = 1719 + 9316 = 11035 kg/m3
// At Pu=0.30 → 11460*0.30 + 10960*0.70 = 3438 + 7672 = 11110 kg/m3
static double moxTheoreticalDensity(double pu) {
    return 11460.0 * pu + 10960.0 * (1.0 - pu);
}

// Melting temperature from FRED.f90 fcp:
//   Tm = 3120 - 388.1*Pu - 30.4*Pu^2  [K]
// At Pu=0   → 3120 K  (pure UO2)
// At Pu=0.15 → 3120 - 58.2 - 0.7 = 3061 K
// At Pu=0.30 → 3120 - 116.4 - 2.7 = 3001 K
static double moxMeltingTemp(double pu) {
    return 3120.0 - 388.1 * pu - 30.4 * pu * pu;
}

MOX::MOX(double pu_content, double rho0)
    : m_pu(pu_content),
      m_tden(moxTheoreticalDensity(pu_content)),
      m_tm(moxMeltingTemp(pu_content))
{
    // Default as-fabricated density: 95% of theoretical density
    m_rho0 = (rho0 > 0.0) ? rho0 : 0.95 * m_tden;
}

// Philipponneau (1992), J. Nucl. Mater. 188 (1992) 194-197.
// FRED.f90 function flamb, zero burnup, stoichiometric fuel (sto=2.0):
//   ac = 1.320*sqrt((2-sto)+0.0093) - 0.091 + 0.0038*bup/9.33
//   At bup=0, sto=2: ac = 1.320*sqrt(0.0093) - 0.091 = 0.036286
//   porosity correction: por = 1 - rho/tden
//   k(T) = [1/(ac + 2.493e-4*T) + 88.4e-12*T^3] * (1-por)/(1+2*por)
double MOX::thermalConductivity(double T) const {
    static const double ac = 1.320 * std::sqrt(0.0093) - 0.091; // 0.036286 at zero burnup, stoichiometric
    double por = 1.0 - m_rho0 / m_tden;
    return (1.0 / (ac + 2.493e-4 * T) + 88.4e-12 * T * T * T)
           * (1.0 - por) / (1.0 + 2.0 * por);
}

// Popov et al., ORNL/TM-2000/351 (2000).
// FRED.f90 function fcp: linear combination of UO2 and PuO2 tabular Cp.
//   cp_MOX(T) = (1-Pu)*cp_UO2(T) + Pu*cp_PuO2(T)
double MOX::heatCapacity(double T) const {
    static const double tTab[] = {
        300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
        1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,
        2500,2600,2700,2800,2900,3000,3100
    };
    static const double cpUO2[] = {
        235.51,265.79,282.14,292.21,299.11,304.24,308.32,311.74,
        314.76,317.59,320.44,323.60,327.41,332.31,338.77,347.30,
        358.40,372.54,390.12,411.46,436.78,466.21,499.80,537.50,
        579.17,624.63,673.63,725.88,781.08
    };
    static const double cpPuO2[] = {
        203.71,225.79,235.77,241.33,244.92,247.47,249.44,251.04,
        252.43,253.70,254.96,256.32,257.93,259.93,262.47,265.67,
        269.57,274.17,279.35,284.96,290.80,296.70,302.47,308.00,
        313.21,318.08,322.59,326.76,330.64
    };
    constexpr int N = 29;

    auto interp = [&](const double* tab) -> double {
        if (T <= tTab[0])   return tab[0];
        if (T >= tTab[N-1]) return tab[N-1];
        int i = 0;
        while (i < N-1 && T > tTab[i+1]) ++i;
        double f = (T - tTab[i]) / (tTab[i+1] - tTab[i]);
        return tab[i] + f * (tab[i+1] - tab[i]);
    };

    return (1.0 - m_pu) * interp(cpUO2) + m_pu * interp(cpPuO2);
}

// MATPRO (NUREG/CR-0479, 1981). FRED.f90 function ftexp.
// eps_UO2(T)  = 1e-5*T - 3e-3 + 4e-2*exp(-5000/T)
// eps_PuO2(T) = 0.9e-5*T - 2.7e-3 + 7e-2*exp(-5072/T)
// eps_MOX(T)  = (1-Pu)*eps_UO2(T) + Pu*eps_PuO2(T)   [absolute strain]
// Returned value = eps_MOX(T) - eps_MOX(293.15) so it is zero at fabrication.
double MOX::thermalExpansionStrain(double T) const {
    auto epsUO2 = [](double t) {
        return 1.0e-5 * t - 3.0e-3 + 4.0e-2 * std::exp(-5000.0 / t);
    };
    auto epsPuO2 = [](double t) {
        return 0.9e-5 * t - 2.7e-3 + 7.0e-2 * std::exp(-5072.0 / t);
    };
    auto epsMOX = [&](double t) {
        return (1.0 - m_pu) * epsUO2(t) + m_pu * epsPuO2(t);
    };
    return epsMOX(T) - epsMOX(T_REF);
}

// MATPRO (NUREG/CR-0479, 1981). FRED.f90 function felmod.
//   tden = 11460*Pu + 10960*(1-Pu)
//   por  = 1 - density/tden
//   E = 2.334e11 * (1 - 2.752*por) * (1 - 1.0915e-4*T) * (1 + 0.15*Pu)  [Pa]
//   Near-melting softening: applied within 150 K of Tm.
//   Note: FRED.f90 original has a sign error in the softening branch
//         (max(1,1-ratio^12) should be min); implemented as in FRED.f90 to match legacy.
double MOX::youngsModulus(double T, double density) const {
    double por = 1.0 - density / m_tden;
    double E   = 2.334e11 * (1.0 - 2.752 * por) * (1.0 - 1.0915e-4 * T)
                 * (1.0 + 0.15 * m_pu);
    double t0 = m_tm - 150.0;
    if (T > t0) {
        double ratio = (T - t0) / (m_tm - t0);
        double fsoft = std::max(0.0, 1.0 - std::pow(ratio, 12.0));
        E *= fsoft;
    }
    return E / 1.0e6; // convert Pa → MPa
}

double MOX::poissonRatio()       const { return 0.276; }
double MOX::referenceDensity()   const { return m_rho0; }
double MOX::theoreticalDensity() const { return m_tden; }
double MOX::meltingTemperature() const { return m_tm; }

} // namespace fred
