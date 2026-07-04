#pragma once
#include "platform/FuelPelletMaterial.hpp"

namespace fred {

// Dummy fuel pellet for benchmark/testing — all properties constant.
//   k   = 10   W/(m·K)
//   rho = 10000 kg/m3
//   cp  = 100  J/(kg·K)
//   CTE = 10e-6 /K  →  ε_th(T) = 1e-5 * (T − 293.15)
//   E   = 100000 MPa  (100 GPa)
//   ν   = 0.3
class DummyFuelPellet : public FuelPelletMaterial {
public:
    double thermalConductivity(double T) const override;
    double heatCapacity(double T) const override;
    double thermalExpansionStrain(double T) const override;
    double youngsModulus(double T, double density) const override;
    double poissonRatio() const override;
    double referenceDensity() const override;
    double theoreticalDensity() const override;
    double meltingTemperature() const override;
};

} // namespace fred
