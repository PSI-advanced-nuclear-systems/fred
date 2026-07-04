#include "FredOxMOX.hpp"
#include <cmath>
#include <algorithm>

// Correlations from legacy FRED.f90: flamb (burnup-dependent), fcp, ftexp, felmod, fpoir.

namespace fred {

static constexpr double T_REF_MOX = 293.15; // K

static double moxTheoreticalDensity(double pu) {
    return 11460.0 * pu + 10960.0 * (1.0 - pu);
}

static double moxMeltingTemp(double pu) {
    return 3120.0 - 388.1 * pu - 30.4 * pu * pu;
}

FredOxMOX::FredOxMOX(double pu_content, double rof0, double sto0)
    : m_pu(pu_content),
      m_rho0(rof0),
      m_tden(moxTheoreticalDensity(pu_content)),
      m_tm(moxMeltingTemp(pu_content)),
      m_sto(sto0)
{}

// Philipponneau (1992) thermal conductivity with burnup degradation term.
// FRED.f90 function flamb:
//   ac = 1.320*sqrt((2-sto)+0.0093) - 0.091 + 0.0038*bup/9.33
//   por = 1 - rho/tden
//   k(T) = [1/(ac + 2.493e-4*T) + 88.4e-12*T^3] * (1-por)/(1+2*por)
double FredOxMOX::thermalConductivity(double T) const {
    const double ac  = 1.320 * std::sqrt((2.0 - m_sto) + 0.0093) - 0.091
                       + 0.0038 * m_bup / 9.33;
    const double por = 1.0 - m_rho0 / m_tden;
    return (1.0 / (ac + 2.493e-4 * T) + 88.4e-12 * T * T * T)
           * (1.0 - por) / (1.0 + 2.0 * por);
}

// Popov et al. ORNL/TM-2000/351: cp = (1-Pu)*cp_UO2(T) + Pu*cp_PuO2(T)
double FredOxMOX::heatCapacity(double T) const {
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

double FredOxMOX::thermalExpansionStrain(double T) const {
    auto epsUO2 = [](double t) {
        return 1.0e-5 * t - 3.0e-3 + 4.0e-2 * std::exp(-5000.0 / t);
    };
    auto epsPuO2 = [](double t) {
        return 0.9e-5 * t - 2.7e-3 + 7.0e-2 * std::exp(-5072.0 / t);
    };
    auto epsMOX = [&](double t) {
        return (1.0 - m_pu) * epsUO2(t) + m_pu * epsPuO2(t);
    };
    return epsMOX(T) - epsMOX(T_REF_MOX);
}

double FredOxMOX::youngsModulus(double T, double density) const {
    double por = 1.0 - density / m_tden;
    double E   = 2.334e11 * (1.0 - 2.752 * por) * (1.0 - 1.0915e-4 * T)
                 * (1.0 + 0.15 * m_pu);
    double t0 = m_tm - 150.0;
    if (T > t0) {
        double ratio = (T - t0) / (m_tm - t0);
        double fsoft = std::max(0.0, 1.0 - std::pow(ratio, 12.0));
        E *= fsoft;
    }
    return E / 1.0e6;
}

double FredOxMOX::poissonRatio()       const { return 0.276; }
double FredOxMOX::referenceDensity()   const { return m_rho0; }
double FredOxMOX::theoreticalDensity() const { return m_tden; }
double FredOxMOX::meltingTemperature() const { return m_tm; }

} // namespace fred
