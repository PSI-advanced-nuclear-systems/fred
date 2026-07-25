#pragma once
#include "apps/fred_m_na/FredMNaCladdingMaterial.hpp"

namespace fred {

// HT-9 ferritic-martensitic steel cladding for FRED-M-Na.
//
// Implements the FredMNaCladdingMaterial interface (CladdingMaterial plus
// FRED-M-Na's irradiation-response virtuals) using correlations from
// Hofmann (1985), Karahan (2007), the SAS4A User's Manual, and the
// Timpano/EPFL FRED_M_OCT24 code.
//
//   Thermal conductivity : Polynomial fit to HT-9 data (Clamb.for)
//   Heat capacity        : ~460 J/(kg·K), approximately constant
//   Thermal expansion    : Karahan (2007) polynomial (Ctexp.for)
//   Young's modulus      : Celmod.for: E = 2.137e11 - 1.0274e8*T  [MPa]
//   Poisson's ratio      : 0.28 (mid-temperature average; T-dependent form
//                          available in PoissonRatio(T) — see HT9.cpp)
//   Meyer hardness       : 3 * sigma_y  (Tabor relation)
//   Reference density    : 7750 kg/m3 (default)
//   Yield stress         : Csigy.for table interpolation [Pa]
//   Creep rate           : Ccreep.for (thermal + irradiation terms) [1/s]
//   Void swelling        : Hofmann (1985), Cswel.for icswel=1 [-]
//   Void swelling (SAS4A): SAS4A User's Manual, Cswel.for icswel=2 [1/s]
//   Wastage (PK/LaTracking): Clanth.for / MFUEL App. 9.3 coefficients
class HT9 : public FredMNaCladdingMaterial {
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

    // FredMNaCladdingMaterial interface
    double yieldStress (double T)                                       const override;
    double creepRate   (double T, double sig_Pa, double qqv, double time_s) const override;
    double voidSwelling(double neuflue, double T)                       const override;
    double voidSwellingSAS4ARate(double T_K, double flux_1e22,
                                  double dconv, double dose_dpa)        const override;
    double wastageDiffusionPK(double T_K)                               const override;
    double wastageSolubilityFuel()                                      const override;
    double wastageSolubilityClad()                                      const override;
    double wastageDiffusionLaTracking(double T_K)                       const override;
    double wastageUptakeCapacity()                                      const override;

    // Burst stress [Pa] — Csigb.for (HT-9 branch), Cfail.for.
    // ssy0 : yield stress [Pa] from yieldStress(T).
    // NOTE: the Fortran Csigb.for contains a bug (tk-tk=0 → constant result).
    //       This implementation uses the corrected temperature-dependent formula:
    //       csigb = ssy0 * (1.1 - 0.1 * tanh((T_Celsius - 200) / 200))
    double burstStress(double T_K, double ssy0)                        const override;

private:
    double m_rho0;
};

} // namespace fred
