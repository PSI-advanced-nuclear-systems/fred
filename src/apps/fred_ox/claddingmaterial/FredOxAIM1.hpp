#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// AIM1 austenitic stainless steel for FRED-OX.
//
// Same elastic/thermal properties as FRED-ROD AIM1 (Luzzi 2013).
// Adds irradiation and thermal creep rate method (FRED.f90 ccreep function).
//
// Thermal creep + irradiation creep model (Luzzi et al., RdS/PAR2013/022):
//   creep_thermal (%/h) = 2.3e14 * exp(-8.46e4 / (R*T)) * sinh(39.72*sigma/(R*T))
//   creep_irrad  (%/h) = 3.2e-24 * E_n [MeV] * phi [n/cm2/s] * sigma [MPa]
//   total creep [1/s]  = (creep_thermal + creep_irrad) / 100 / 3600
// with R=1.986 cal/(mol·K), E_n=2 MeV, phi=1e18 n/(cm2·s).
class FredOxAIM1 : public CladdingMaterial {
public:
    explicit FredOxAIM1(double rho0 = 7900.0) : m_rho0(rho0) {}

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T) const override;
    double poissonRatio() const override;
    double meyerHardness(double T) const override;
    double referenceDensity() const override;

    // Effective creep strain rate [1/s] at temperature T [K] and hoop stress sigma [MPa].
    double creepRate(double T, double sigma) const override;

private:
    double m_rho0;
};

} // namespace fred
