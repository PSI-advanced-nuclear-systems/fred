#pragma once
#include "platform/FuelPelletMaterial.hpp"

namespace fred {

// MOX fuel pellet for FRED-OX — includes burnup-dependent thermal conductivity.
//
// Extends the FRED-ROD MOX by adding local burnup and stoichiometry state
// that modify the Philipponneau (1992) thermal conductivity formula.
// All other properties (cp, thermal expansion, Young's modulus) are the same
// as the fresh-fuel MOX correlation.
//
// The solver calls setBurnup() and setStoichiometry() before each residual
// evaluation so that thermalConductivity() returns the correct value.
//
// Thermal conductivity (Philipponneau 1992 + burnup degradation):
//   ac = 1.320 * sqrt((2-sto) + 0.0093) - 0.091 + 0.0038*bup/9.33
//   k(T) = [1/(ac + 2.493e-4*T) + 88.4e-12*T^3] * (1-por)/(1+2*por)
//
// Heat capacity:  Popov et al. ORNL/TM-2000/351
// Thermal expansion: MATPRO
// Young's modulus: MATPRO with porosity correction
class FredOxMOX : public FuelPelletMaterial {
public:
    // pu_content : Pu mole fraction [-], e.g. 0.18
    // rof0       : as-fabricated density [kg/m3]
    // sto0       : initial stoichiometry (2.0 = stoichiometric; typically 1.95-1.99)
    explicit FredOxMOX(double pu_content, double rof0, double sto0 = 1.97);

    // Update irradiation state (called by solver before each residual evaluation)
    void setBurnup(double bup_MWdkgU)  { m_bup = bup_MWdkgU; }
    void setStoichiometry(double sto)  { m_sto = sto; }
    double burnup()      const { return m_bup; }
    double puContent()   const { return m_pu; }
    double stoichiometry() const { return m_sto; }

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T, double density) const override;
    double poissonRatio() const override;
    double referenceDensity() const override;
    double theoreticalDensity() const override;
    double meltingTemperature() const override;

private:
    double m_pu;    // Pu mole fraction [-]
    double m_rho0;  // as-fabricated density [kg/m3]
    double m_tden;  // theoretical density [kg/m3]
    double m_tm;    // melting temperature [K]
    double m_bup = 0.0; // current burnup [MWd/kgU]
    double m_sto;       // current stoichiometry
};

} // namespace fred
