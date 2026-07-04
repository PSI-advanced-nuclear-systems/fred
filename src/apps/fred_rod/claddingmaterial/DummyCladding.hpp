#pragma once
#include "platform/CladdingMaterial.hpp"

namespace fred {

// Dummy cladding for benchmark/testing — all properties constant.
//   k   = 10   W/(m·K)
//   rho = 10000 kg/m3
//   cp  = 100  J/(kg·K)
//   CTE = 10e-6 /K  →  ε_th(T) = 1e-5 * (T − 293.15)
//   E   = 100000 MPa  (100 GPa)
//   ν   = 0.3
//   H_M = 1e9 Pa  (Meyer hardness)
class DummyCladding : public CladdingMaterial {
public:
    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T) const override;
    double poissonRatio() const override;
    double meyerHardness(double T) const override;
    double referenceDensity() const override;
};

} // namespace fred
