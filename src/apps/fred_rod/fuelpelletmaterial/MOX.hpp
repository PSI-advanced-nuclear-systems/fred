#pragma once
#include "platform/FuelPelletMaterial.hpp"

namespace fred {

// Mixed-oxide (MOX) fuel pellet — FRED-ROD implementation.
//
// Valid for fresh (unirradiated) fuel at as-fabricated density.
// All correlations taken directly from legacy FRED.f90 (flamb, fcp, ftexp,
// felmod, fpoir functions), which implement:
//
//   Thermal conductivity : Y. Philipponneau, "Thermal conductivity of
//     (U,Pu)O2-x mixed oxide fuel", J. Nucl. Mater. 188 (1992) 194-197.
//     At zero burnup and stoichiometric fuel (sto = 2.0):
//       ac = 1.320*sqrt(0.0093) - 0.091 = 0.03629
//       k(T) = [1/(ac + 2.493e-4*T) + 88.4e-12*T^3] * (1-por)/(1+2*por)
//     Porosity uses Pu-dependent theoretical density: rho_t = 11460*Pu + 10960*(1-Pu).
//
//   Heat capacity  : S.G. Popov et al., "Thermophysical properties of MOX and
//     UO2 fuels including the effects of irradiation", ORNL/TM-2000/351 (2000).
//     Tabular data for UO2 and PuO2 cp; linear combination: cp = (1-Pu)*cp_UO2 + Pu*cp_PuO2.
//
//   Thermal expansion : MATPRO (MATPRO-09, 1981; NUREG/CR-0479).
//     eps_UO2(T)  = 1e-5*T - 3e-3 + 4e-2*exp(-5000/T)
//     eps_PuO2(T) = 0.9e-5*T - 2.7e-3 + 7e-2*exp(-5072/T)
//     eps_MOX(T)  = (1-Pu)*eps_UO2(T) + Pu*eps_PuO2(T)   (strain relative to 293.15 K)
//
//   Young's modulus : MATPRO with porosity correction and Pu enhancement factor (1+0.15*Pu).
//     E(T,rho) = 2.334e11 * (1 - 2.752*por) * (1 - 1.0915e-4*T) * (1 + 0.15*Pu)
//     Near-melting softening applied within 150 K of Tm.
//
//   Poisson's ratio : 0.276 (MATPRO).
//
//   Melting temperature : Tm = 3120 - 388.1*Pu - 30.4*Pu^2  [K]
//     (MATPRO correlation for MOX liquidus).
//
// Realistic Pu content values for fast reactor fuel pins:
//   Low  (MOX_low)  : pu_content = 0.15  (15 mol%)  — BN-800 / EFR inner pins
//   High (MOX_high) : pu_content = 0.30  (30 mol%)  — BN-800 driver fuel / outer core pins
//
// Example usage:
//   fred::MOX mox_low (0.15);   // pu_content = 0.15
//   fred::MOX mox_high(0.30);   // pu_content = 0.30
class MOX : public FuelPelletMaterial {
public:
    // pu_content : plutonium mole fraction [-], e.g. 0.15 or 0.30
    // rho0       : as-fabricated density [kg/m3]; if <=0, defaults to 95% of
    //              the Pu-dependent theoretical density.
    explicit MOX(double pu_content, double rho0 = -1.0);

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T, double density) const override;
    double poissonRatio() const override;
    double referenceDensity() const override;
    double theoreticalDensity() const override;
    double meltingTemperature() const override;

private:
    double m_pu;    // plutonium mole fraction [-]
    double m_rho0;  // as-fabricated density [kg/m3]
    double m_tden;  // theoretical density [kg/m3]
    double m_tm;    // melting temperature [K]
};

} // namespace fred
