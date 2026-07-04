#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// T91 (9Cr-1Mo-V) ferritic-martensitic steel cladding — FRED-ROD implementation.
//
// No irradiation effects (creep, swelling, embrittlement) — these would be added
// in a FRED-OX-style derived class.
//
// All correlations taken from legacy FRED.f90 (clamb, ccp, ctexp, celmod, cpoir
// functions for cmat='t91'), which implement:
//
//   Thermal conductivity : k [W/(m·K)], Tc in °C
//     k = 23.71 + 0.01718*Tc - 1.45e-5*Tc^2
//     Source: FRED.f90 clamb; referenced against OECD/NEA data for T91.
//     (Consistent with Leibowitz & Blomquist, Int. J. Thermophysics 9(5):873-883, 1988,
//      and Natesan et al., NUREG/CR-6612 data for 9Cr-1Mo-V steels.)
//
//   Heat capacity : cp [J/(kg·K)], T in K
//     cp = 431.0 + 0.177*T + 8.72e-5*T^2
//     Source: FRED.f90 ccp (same formula as AIM1); attributed to
//     L. Luzzi et al., "Modeling and Analysis of Nuclear Fuel Pin Behavior for
//     Innovative Lead Cooled FBR", Report RdS/PAR2013/022.
//
//   Thermal expansion strain : absolute strain [-], T in K
//     eps(T) = -0.2177e-2 + 6.735e-6*T + 5.12e-9*T^2 - 2.248e-12*T^3 + 3.933e-16*T^4
//     Strain returned = eps(T) - eps(293.15).
//     Source: FRED.f90 ctexp (same polynomial as AIM1); Luzzi et al. RdS/PAR2013/022.
//
//   Young's modulus : E [MPa], Tc in °C
//     Tc <= 500°C : E = (2.073e11 - 6.458e7*Tc)  [Pa]
//     Tc >  500°C : E = (2.95e11  - 2.40e8*Tc)   [Pa]
//     Source: FRED.f90 celmod for cmat='t91'.
//
//   Poisson's ratio : 0.28
//     Source: FRED.f90 cpoir for cmat='t91'.
//
//   Meyer hardness : SAS4A correlation (same as AIM1), used for gap contact conductance.
//     Source: FRED.f90 gaphtc, SAS4A/SASSYS-1 hgap subroutine.
class T91 : public CladdingMaterial {
public:
    // rho0: as-fabricated density [kg/m3]. T91 nominal: 7750 kg/m3.
    explicit T91(double rho0 = 7750.0);

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
