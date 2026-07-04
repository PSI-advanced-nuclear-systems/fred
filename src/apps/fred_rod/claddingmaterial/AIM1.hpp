#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// AIM1 austenitic stainless steel cladding.
//
// FRED-OX note: FredOxCladdingMaterial will inherit from this class to add
// irradiation creep, void swelling, and embrittlement correlations.
//
// Source: L. Luzzi et al., "Modeling and Analysis of Nuclear Fuel Pin Behavior
//         for Innovative Lead Cooled FBR", Report RdS/PAR2013/022.
//
//   Thermal conductivity : k = 13.95 + 0.01163 * Tc  [W/(m·K)],  Tc in °C
//   Heat capacity        : cp = 431.0 + 0.177*T + 8.72e-5*T^2  [J/(kg·K)]
//   Thermal expansion    : Luzzi polynomial in T [K]
//   Young's modulus      : E = 2.027e11 - 8.167e7*Tc  [Pa]
//   Poisson's ratio      : 0.289
//   Meyer hardness       : SAS4A correlation (for gap contact conductance)
class AIM1 : public CladdingMaterial {
public:
    explicit AIM1(double rho0 = 7900.0);

    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T) const override;
    double poissonRatio() const override;
    double meyerHardness(double T) const override;
    double referenceDensity() const override;

private:
    double m_rho0;
};

} // namespace fred
