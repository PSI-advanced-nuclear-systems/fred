#pragma once
#include "platform/FuelPelletMaterial.hpp"

namespace fred {

// Uranium dioxide (UO2) fuel pellet — FRED-ROD implementation.
//
// Valid for fresh (unirradiated) fuel at as-fabricated density.
// Thermal conductivity uses the Philipponneau (1992) correlation at
// zero burnup and as-fabricated porosity.
//
// FRED-OX note: FredOxFuelPelletMaterial inherits from this class and
// overrides thermalConductivity to accept current density and burnup,
// enabling the full irradiation-dependent Philipponneau correlation.
//
// References:
//   Thermal conductivity : Philipponneau (1992), J. Nucl. Mat. 188, 194-197
//   Heat capacity        : Popov et al. ORNL/TM-2000/351
//   Thermal expansion    : MATPRO
//   Young's modulus      : MATPRO with porosity and near-melting softening
//   Poisson's ratio      : MATPRO, 0.276
class UO2 : public FuelPelletMaterial {
public:
    // rho0: as-fabricated density [kg/m3].  Theoretical density of UO2 is 10960 kg/m3.
    explicit UO2(double rho0 = 10400.0);

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T, double density) const override;
    double poissonRatio() const override;
    double referenceDensity() const override;
    double theoreticalDensity() const override;
    double meltingTemperature() const override;

private:
    double m_rho0;
};

} // namespace fred
