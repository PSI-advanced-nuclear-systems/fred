#include "FredOxGapMaterial.hpp"
#include "apps/fred_ox/FredOxIrradiationPhysics.hpp"
#include <cmath>
#include <algorithm>

namespace fred {

// ---------------------------------------------------------------------------
// Gas mixture conductivity — moved from FredOxIrradiationPhysics (material
// property of the gap, not a general irradiation physics function).
// ---------------------------------------------------------------------------
double FredOxGapMaterial::gasMixtureConductivity(double mu0_he, double mu0_ar,
                                                  double fgrel, double T_avg)
{
    const double k_He = 2.639e-3 * std::pow(T_avg, 0.7085);
    const double k_Ar = 2.986e-4 * std::pow(T_avg, 0.7224);
    const double k_Kr = 8.247e-5 * std::pow(T_avg, 0.8363);
    const double k_Xe = 4.351e-5 * std::pow(T_avg, 0.8618);

    // Mole fractions: fill gas + released fission gas
    // FG composition from legacy FRED gaphtc: 88.46% Xe, 7.69% Kr, 3.85% He
    double x_He = mu0_he + 0.0385 * fgrel;
    double x_Ar = mu0_ar;
    const double x_Kr = 0.0769 * fgrel;
    const double x_Xe = 0.8846 * fgrel;

    const double x_total = x_He + x_Ar + x_Kr + x_Xe;
    if (x_total <= 0.0) return k_He; // pure helium fallback

    // Geometric mean: k_mix = prod(k_i^(x_i/x_total))
    const double ln_k = (x_He * std::log(k_He)
                       + x_Ar * std::log(std::max(k_Ar, 1e-30))
                       + x_Kr * std::log(k_Kr)
                       + x_Xe * std::log(k_Xe)) / x_total;
    return std::exp(ln_k);
}

// ---------------------------------------------------------------------------
// GapMaterial interface
// ---------------------------------------------------------------------------
double FredOxGapMaterial::gapConductivity(double T) const {
    return gasMixtureConductivity(m_mu0, 0.0, m_fgrel, T);
}

double FredOxGapMaterial::relocationFraction() const {
    if (!m_reloc || m_burnup <= 0.0) return 0.0;
    return fuelRelocFraction(m_ql, m_burnup);
}


} // namespace fred
