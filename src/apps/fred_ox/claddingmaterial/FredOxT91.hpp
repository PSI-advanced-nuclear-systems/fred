#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// T91 ferritic-martensitic steel for FRED-OX.
//
// Same elastic/thermal properties as FRED-ROD T91.
// Creep rate is zero (legacy FRED.f90 ccreep returns 0 for T91;
// no irradiation creep model is implemented for T91 in the reference code).
class FredOxT91 : public CladdingMaterial {
public:
    explicit FredOxT91(double rho0 = 7750.0) : m_rho0(rho0) {}

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T) const override;
    double poissonRatio() const override;
    double meyerHardness(double T) const override;
    double referenceDensity() const override;

    // Creep rate for T91 is zero in legacy FRED (no irradiation/thermal creep model).
    double creepRate(double /*T*/, double /*sigma*/) const override { return 0.0; }

private:
    double m_rho0;
};

} // namespace fred
