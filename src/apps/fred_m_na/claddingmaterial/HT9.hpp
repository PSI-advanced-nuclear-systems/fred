#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// HT-9 ferritic-martensitic steel cladding for FRED-M-Na.
//
// Implements the CladdingMaterial interface using correlations from
// Hofmann (1985), Karahan (2007), and the Timpano/EPFL FRED_M_OCT24 code.
//
//   Thermal conductivity : Polynomial fit to HT-9 data (Clamb.for)
//   Heat capacity        : ~460 J/(kg·K), approximately constant
//   Thermal expansion    : Karahan (2007) polynomial (Ctexp.for)
//   Young's modulus      : Celmod.for: E = 2.137e11 - 1.0274e8*T  [MPa]
//   Poisson's ratio      : 0.28 (mid-temperature average; T-dependent form
//                          available in HT9::poissonRatio(T) — see HT9.cpp)
//   Meyer hardness       : 3 * sigma_y  (Tabor relation)
//   Reference density    : 7750 kg/m3 (default)
//
// Extended methods beyond CladdingMaterial (accessible via HT9&):
//   yieldStress(T)                       — Csigy.for table interpolation [Pa]
//   creepRate(T, sig_Pa, qqv, time_s)    — Ccreep.for [1/s]
//   voidSwelling(neuflue, T)             — Cswel.for model 1 [-]
class HT9 : public CladdingMaterial {
public:
    explicit HT9(double rho0 = 7750.0);

    // CladdingMaterial interface
    double thermalConductivity  (double T) const override;
    double heatCapacity         (double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus        (double T) const override;
    double poissonRatio         ()         const override;
    double meyerHardness        (double T) const override;
    double referenceDensity     ()         const override;

    // HT9-specific extended correlations
    double yieldStress (double T)                                       const;
    double creepRate   (double T, double sig_Pa, double qqv, double time_s) const;
    double voidSwelling(double neuflue, double T)                       const;

    // Burst stress [Pa] — Csigb.for (HT-9 branch), Cfail.for.
    // ssy0 : yield stress [Pa] from yieldStress(T).
    // NOTE: the Fortran Csigb.for contains a bug (tk-tk=0 → constant result).
    //       This implementation uses the corrected temperature-dependent formula:
    //       csigb = ssy0 * (1.1 - 0.1 * tanh((T_Celsius - 200) / 200))
    double burstStress(double T_K, double ssy0)                        const;

private:
    double m_rho0;
};

} // namespace fred
