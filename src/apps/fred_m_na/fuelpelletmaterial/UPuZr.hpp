#pragma once
#include "platform/FuelPelletMaterial.hpp"

namespace fred {

// Irradiation correction model for U-Pu-Zr thermal conductivity.
// Corresponds to the integer flag `f` in Flamb.for (Timpano/EPFL FRED_M_OCT24).
enum class ConductivityModel {
    DetailedNaSodium = 1, // f=1  Na-infiltration Maxwell-Eucken (Karahan 2009 + MFUEL)
    EmpiricalBurnup  = 2, // f=2  piecewise empirical burnup degradation
    EsfrSimple       = 3  // f=3  ESFR-SIMPLE sigmoid fit (requires only burnup)
};

// U-Pu-Zr metallic fuel pellet for FRED-M-Na.
//
// Implements the FuelPelletMaterial interface using correlations from
// Karahan (2009) and the Timpano/EPFL FRED_M_OCT24 Fortran code.
//
// Fresh-fuel properties (no irradiation state in the base class):
//   Thermal conductivity   : fresh-fuel Aydin/Karahan formula (no burnup correction)
//                            Use UPuZr::thermalConductivityIrradiated() for
//                            the irradiation-corrected version with Na infiltration.
//   Heat capacity          : Karahan (2009), Hales (2016) correction
//   Thermal expansion      : Karahan (2009)
//   Young's modulus        : Karahan (2009) + porosity correction
//   Poisson's ratio        : Karahan (2009) + porosity correction
//   Reference density      : as-fabricated [kg/m3]
//   Theoretical density    : calculated from composition
//   Melting temperature    : approximate (Tm ~ 1450 K for alpha, 1378 K for full U-Pu-Zr)
//
// Composition:
//   pu_weight_frac : Pu weight fraction [-]  (typical: 0.10 – 0.22 for SFR pins)
//   zr_weight_frac : Zr weight fraction [-]  (typical: 0.06 – 0.10)
//   U is the remainder.
class UPuZr : public FuelPelletMaterial {
public:
    // pu_weight_frac : Pu weight fraction [-]
    // zr_weight_frac : Zr weight fraction [-]
    // rho0           : as-fabricated density [kg/m3]; if <=0 defaults to 75% of theoretical
    explicit UPuZr(double pu_weight_frac, double zr_weight_frac, double rho0 = -1.0);

    // Conductivity model selector (default: DetailedNaSodium).
    void              setConductivityModel(ConductivityModel m) { m_conductivity_model = m; }
    ConductivityModel conductivityModel()                  const { return m_conductivity_model; }

    // --- FuelPelletMaterial interface ---

    // Fresh-fuel thermal conductivity (no irradiation correction) [W/(m·K)]
    double thermalConductivity(double T) const override;

    // Irradiation-corrected thermal conductivity [W/(m·K)].
    // Dispatches to the model selected by setConductivityModel():
    //   DetailedNaSodium — Maxwell-Eucken with Na infiltration (needs poros_tot/gas, psod)
    //   EmpiricalBurnup  — piecewise burnup degradation (needs bup_FIMA)
    //   EsfrSimple       — ESFR-SIMPLE sigmoid fit (needs bup_FIMA)
    double thermalConductivityIrradiated(double T_K,
                                          double bup_FIMA,
                                          double poros_tot,
                                          double poros_gas,
                                          double psod) const;

    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T, double density) const override;
    double poissonRatio() const override;
    double referenceDensity() const override;
    double theoreticalDensity() const override;
    // Solidus temperature for the as-fabricated composition [K].
    // Uses the 25-point Zr-dependent table from Mmelt.for (Pu=13.25 wt% fixed).
    // Nearest-neighbour lookup on Zr weight-percent (5–29 wt%).
    double meltingTemperature() const override;

    // Solidus temperature [K] for an arbitrary Zr weight fraction (0–1).
    // Same 25-point table lookup; use when Zr content varies by node.
    double solidusTemperature(double zr_wf) const;

    // Accessors for fuel composition
    double puContent() const { return m_pu; }
    double zrContent() const { return m_zr; }

    // Per-node overloads that accept explicit local composition (pu_wf, zr_wf)
    // rather than the as-fabricated global m_pu/m_zr.  Used by
    // FredMNaSolver::afterAcceptedStep to compute k_irr_factor from the
    // post-redistribution Zr/Pu profile instead of the nominal composition.
    double thermalConductivityLocal(double T, double pu_wf, double zr_wf) const;
    double thermalConductivityIrradiatedLocal(double T_K, double pu_wf, double zr_wf,
                                               double bup_FIMA,
                                               double poros_tot, double poros_gas,
                                               double psod) const;

    // UPuZr-specific extended correlations
    // Full temperature- and density-dependent Poisson's ratio [-]
    double poissonRatioFull(double T_K, double density) const;
    // Effective creep strain rate [1/s] — Karahan (2009) / Hofman / Gruber
    //   sig_eff    : von Mises effective stress [Pa]
    //   sgh,sgz,sgr: hoop, axial, radial stress components [Pa]
    //   gasp_Pa    : internal gas pressure [Pa]
    //   dotfissden : fission rate density [fiss/m3/s]
    //   swopen     : open-pore swelling fraction (GRSIS; 0 for simplified)
    double creepRate(double T_K, double sig_eff,
                     double sgh, double sgz, double sgr,
                     double gasp_Pa, double dotfissden,
                     double swopen = 0.0) const;

private:
    double m_pu;   // Pu weight fraction [-]
    double m_zr;   // Zr weight fraction [-]
    double m_rho0; // as-fabricated density [kg/m3]
    double m_tden; // theoretical density [kg/m3]
    ConductivityModel m_conductivity_model = ConductivityModel::DetailedNaSodium;
};

} // namespace fred
